# nearscaff

参考基因组引导的基因组挂载工具。`nearscaff` 利用近缘参考基因组的
蛋白序列和核苷酸共线性，把 contig 级别的 query 组装锚定到参考基因组上，
并对 contig 进行排序和定向，输出染色体级别的 scaffold。

## 工作原理

- **Stage 0 —— 蛋白锚定与 Block Tree 构建。**
  用 [miniprot](https://github.com/lh3/miniprot) 把参考蛋白（直接提供，
  或从参考 GFF3 + 基因组中提取）比对到 query 组装上。锚点经 C-score
  过滤后，再用 minimap2 把 query contig 比对到参考基因组以获得真实的
  基因组坐标，随后把 contig 级别的共线性区块组装成 *Block Tree*
  （亚基因组 → 染色体 → 局部簇的层级结构）。输出：`block_tree.json`。
- **Stage 1 —— 挂载（scaffolding）。** 以 Block Tree 中的蛋白-共线性
  边为种子构建 scaffold 图，再通过多轮约束性 minimap2 核苷酸比对
  （`asm5` → `asm10` → `asm20`）逐步把未挂载的 contig 拉入 scaffold。
  最终 contig 顺序由最大权匹配确定（小图用 Edmonds Blossom 精确算法，
  边数超过 20,000 时切换为贪心启发式）。scaffold 按染色体合并后写出
  AGP，并转换为 FASTA。

## 依赖

Python ≥ 3.10，以及 Python 包：

- `networkx`

PATH 中需要的外部工具：

- [minimap2](https://github.com/lh3/minimap2)
- [miniprot](https://github.com/lh3/miniprot)
- `samtools`（用于提取参考基因组区域）

## 安装

```bash
pip install .
```

## 用法

完整流程（Stage 0 + Stage 1）：

```bash
nearscaff run -r REF.fa -p PROTEINS.pep -q QUERY.fa -o OUT -t 16 --preset near
```

也可以用 `-g REF.gff3` 代替 `-p`，从参考 GFF3 注释中提取蛋白。常用选项：

- `--preset {close,near,mid,far}` —— 分化程度预设；决定 miniprot 剪接
  模型和 minimap2 预设（`near` → `asm10`，以此类推）
- `--cscore FLOAT` —— 锚点 C-score 阈值（默认 0.7）
- `--min-cluster-size INT` —— 每个共线性区块的最少锚点数（默认 4）
- `--overlap-threshold FLOAT` —— 区块重叠阈值（默认 0.5）
- `--no-best-buddy` —— 关闭图求解中的 best-buddy 权重缩放
- `--nucleotide-passes asm5 asm10 asm20` —— 核苷酸扩展轮次
- `--keep-intermediate` —— 保留中间文件（锚点 TSV、区块文件）

只跑 Stage 1（基于已有的 Block Tree）：

```bash
nearscaff scaffold -b OUT/block_tree.json -r REF.fa -q QUERY.fa -o OUT2 -t 16
```

Stage 1 额外选项：`--margin`（比对区域外延长度，默认 50000）和
`--unknown-gap-size`（默认 100）。

## 输出

输出目录中的文件：

- `block_tree.json` —— Stage 0 生成的层级 Block Tree
- `nearscaff.agp` —— 最终 scaffold，AGP v2.1 格式
- `nearscaff.scaffolds.fa` —— 最终 scaffold 序列。未能挂载的 query
  contig 会以原始（完整）FASTA 头追加在文件末尾
- `nearscaff_tiered.paf` —— 已挂载 contig 的类 PAF 格式文件，带置信度
  层级标签（`nc:Z:protein|asm5|asm10|asm20`）
- `nearscaff.log` —— 流程日志
- `ref_proteins.faa`、`query.mpi` —— 蛋白集合与 miniprot 索引
- 使用 `--keep-intermediate` 时：`gene_anchors.tsv`、`synteny.blocks`

## 实际案例

两个 contig 级别的植物组装，使用同一套近缘染色体级别参考基因组
（14 条染色体）和参考蛋白集进行挂载，`--preset near`。

案例 1 —— 高度碎片化的组装：

| | 挂载前 | 挂载后 |
|---|---|---|
| 序列数 | 551,915 条 contig | 14 条染色体级 scaffold |
| 总大小 | 564.0 Mb | 挂载 431.8 Mb（76.6%） |
| N50 | 0.019 Mb | 30.8 Mb（scaffold） |
| BUSCO (embryophyta_odb10) | C:83.0% [S:80.7%, D:2.3%], F:13.3%, M:3.8% | C:85.0% [S:82.9%, D:2.1%], F:10.5%, M:4.5% |

案例 2 —— 长读长 contig 组装：

| | 挂载前 | 挂载后 |
|---|---|---|
| 序列数 | 17,594 条 contig | 14 条染色体级 scaffold |
| 总大小 | 456.0 Mb | 挂载 437.7 Mb（96.0%） |
| N50 | 2.95 Mb | 30.9 Mb（scaffold） |
| BUSCO (embryophyta_odb10) | C:97.7% [S:95.4%, D:2.3%], F:2.1%, M:0.2% | C:97.4% [S:95.4%, D:2.0%], F:1.9%, M:0.7% |

两个案例中 contig 都被挂载为与参考染色体一一对应的 scaffold，挂载后
完整 BUSCO 数量保持或有所提升。

## 开发

运行测试（当 miniprot/minimap2/samtools 不在 PATH 中时，集成测试会
自动跳过）：

```bash
python -m pytest tests -q
```
