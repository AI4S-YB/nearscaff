# nearscaff

参考基因组引导的基因组挂载工具。`nearscaff` 利用近缘参考基因组的
蛋白序列和核苷酸共线性，把 contig 级别的 query 组装锚定到参考基因组上，
并对 contig 进行排序和定向，输出染色体级别的 scaffold。
另附独立的亚基因组分型模块 `kmer-phase`（见下文）。

## 工作原理

- **Stage 0 —— 蛋白锚定与 Block Tree 构建。**
  用 [miniprot](https://github.com/lh3/miniprot) 把参考蛋白（直接提供，
  或从参考 GFF3 + 基因组中提取）比对到 query 组装上。锚点经 C-score
  过滤后，再用 minimap2 把 query contig 比对到参考基因组以获得真实的
  基因组坐标，随后把 contig 级别的共线性区块组装成 *Block Tree*
  （亚基因组 → 染色体 → 局部簇的层级结构）。输出：`block_tree.json`。
- **Stage 1 —— 挂载（scaffolding）。** 以 Block Tree 中的蛋白-共线性
  边为种子构建 scaffold 图，再通过多轮 minimap2 核苷酸比对
  （默认 `asm5` → `asm20` 两轮；中间的 `asm10` 轮是冗余的，因为更宽松
  的 `asm20` 比对结果是它的超集）逐步把未挂载的 contig 拉入 scaffold。
  每轮扩展只对可复用的、按预设区分的参考索引（`.mmi`，参考文件变更时
  自动重建）做一次全参考比对，并保留至多 `-N 5` 个次级比对——主比对
  离所有 scaffold 都远的 contig 也可以通过次级比对挂到 scaffold 末端。
  contig 邻接关系由最大权匹配确定（小图用 Edmonds Blossom 精确算法，
  边数超过 20,000 时切换为贪心启发式）。scaffold 按染色体合并后，
  最终精修直接复用 Stage 0 和扩展轮缓存的逐 contig 比对（只对无缓存的
  contig 重新比对，且用 `--secondary=no`），按参考基因组中点坐标
  排序、按比对链向定向（仅由低 mapq 比对支持的方向标记为 `?`），
  最后写出 AGP 并转换为 FASTA。

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
- `--keep-unplaced [MINLEN]` —— 把从未进入 scaffold 的 query contig
  以原名 singleton 行追加进 AGP（默认关闭；带 MINLEN 时只保留
  ≥MINLEN bp 的）。高度碎片化的 query 若不开启，相当比例的组装
  序列会从输出中丢失，悄悄拉低下游完整性指标
- `--nucleotide-passes asm5 asm20` —— 核苷酸扩展轮次（默认
  `asm5 asm20`；传入 `asm5 asm10 asm20` 可还原 0.3.1 之前的行为）
- `--secondary-alignments INT` —— 扩展轮保留的 minimap2 `-N` 次级
  比对数（默认 5；重复序列多的基因组可调低换速度，代价是少挽回
  少量 contig）
- `--no-reuse-index` —— 不用可复用的 minimap2 参考索引，回退到旧的
  逐区域比对路径
- `--keep-intermediate` —— 保留中间文件（锚点 TSV、区块文件）

只跑 Stage 1（基于已有的 Block Tree）：

```bash
nearscaff scaffold -b OUT/block_tree.json -r REF.fa -q QUERY.fa -o OUT2 -t 16
```

> **关于可复现性：** scaffold 拓扑存在一定的运行间随机性——权重相同
> 的边在图求解时按 hash 顺序破局，因此 scaffold 跨度（和 N50）在不同
> 运行间可能有所波动，但挂载的 contig 集合基本不变。如果某次结果不
> 理想，多跑一两次通常能得到更优的解。

Stage 1 额外选项：`--margin`（比对区域外延长度，默认 50000）、
`--unknown-gap-size`（默认 100）、`--keep-unplaced [MINLEN]`（同上），
以及上文提到的 `--nucleotide-passes`、`--secondary-alignments`、
`--no-reuse-index`。

## 用长读长填充 gap（gapfill）

`nearscaff gapfill` 用长读长（HiFi/ONT）选择性地关闭 scaffold 中的
gap。只有两翼都是高置信 contig 的 gap 才会被填充——置信层级取自
`nearscaff_tiered.paf`，默认为 `protein`/`asm5`；触及低置信
（`asm20`）挂载的 gap 保持为 N。填充只替换 N 串、绝不改变
scaffold 拓扑。

工作原理：为每个可填 gap 的两翼各 ~1.5 kb flank 建一个迷你参考，
用 minimap2 把长读长比对上去，两种闭合方式——横跨两翼的 reads
（取跨越段的中位长度 consensus），以及两翼 clip 尾巴互相 overlap
（dovetail 重叠，且重叠段必须唯一比对到 scaffold 上的单一基因座，
以此过滤重复序列造成的假闭合）。

```bash
nearscaff gapfill -a OUT/nearscaff.agp -q QUERY.fa \
    -r hifi.fastq.gz -o OUT_gapfill -t 16
```

选项：`--lr-preset {map-hifi,map-ont}`（默认 `map-hifi`）；
`--lr-min-span INT` —— 每个 span 闭合的最少横跨 reads（默认 1）；
`--max-depth-factor FLOAT` —— 跨越深度超过中位数 N 倍的 gap 拒绝
闭合（重复塌缩；默认 3.0，0 关闭）；`--no-overlap-closure` —— 只做
span 闭合；`--fill-tiers` —— 可填的两翼层级（默认 `protein,asm5`）。
需要 PATH 中有 `minimap2`。

输出：`nearscaff.gapfill.agp` 与 `nearscaff.gapfill.scaffolds.fa`，
以及填充报告（可填 / 闭合数、填补碱基数，区分 span 与 overlap）。

### 短读长模式（`--method sr`）

也支持二代双端数据（`-1/-2`）。三种机制：单 read 横跨闭合（微小
gap 的精确序列）、双端尾巴 overlap 闭合（拼接后的配对序列，不做
唯一性过滤）、PE 双端 gap 缩放（AGP gap 长度替换为基于插入片段的
估计值，类型 `N`、证据 `paired-ends`；`--sr-insert` 设定名义插入
片段，默认 500）。

> **注意——假阳性率高。** contig 末端是串联重复，短读长证据在
> 那里被多拷贝假象主导。在真实碎片化组装上，sr 模式的序列闭合率
> 约为可填 gap 的 0.1~0.6%（另有数百个 gap 获得尺寸估计）。请为
> 闭合率而用它，不要为正确率；有长读长时优先长读长。

注：短读长盲搜方案（abyss-sealer Bloom 游走、clip 延伸）经过实测
已被移除——真实数据上 0 闭合。锚定在串联重复上的 gap 用任何读长
都无法可靠填充；可闭合的主要集中在真实 gap 较小或连接处为唯一
序列的情形。

## 用转录本填充基因区 gap（rna-fill）

`nearscaff rna-fill` 闭合对注释最重要的 gap：基因区内的那些。
scaffold 化后的基因组 BUSCO 可以很高，但基因区仍然布满 N 串，导致
下游注释（miniprot 等）只能得到破碎的基因模型。rna-fill 把转录本
序列比对回 scaffold，填补转录本能够跨越的 gap。

输入转录本可以是 Iso-Seq、ONT cDNA 或组装好的转录本
（FASTA/FASTQ）。短读长 RNA-seq 必须先组装（如 Trinity）——单条
短读长无法跨越 gap。

```bash
nearscaff rna-fill -a OUT/nearscaff.agp -q QUERY.fa \
    -T transcripts.fa -o OUT_rnafill -t 16
```

**定向招募 reads 再做组装。** 只有二代 RNA-seq 时，可以只组装
"有用的" reads：`--reads` + `--proteins` 会构建组合诱饵（gap flank
——含 `--internal-n` 炸开的组件内 N run——加上**残缺基因位点**：
miniprot 预测出的翻译含 X 的 mRNA，即 CDS 穿过 N 串的位点）并提取
招募 reads 及其 mate：

```bash
nearscaff rna-fill -a OUT/nearscaff.agp -q QUERY.fa \
    --reads RNA_R1.fq.gz RNA_R2.fq.gz --proteins REF.pep.fa \
    --internal-n 20 -o OUT_rnafill -t 16
# 用任意组装器（如 Trinity）组装 OUT_rnafill/rnafill.recruit_*.fastq，
# 然后填充：
nearscaff rna-fill -a OUT/nearscaff.agp -q QUERY.fa \
    -T assembled.fa -o OUT_rnafill -t 16 \
    --internal-n 20 --fill-mode exon-only --one-sided
```

（`--internal-n` 建议搭配 `--fill-mode exon-only`：组件内 N 串多数
盖住内含子，cdna 模式会把内含子压没，扭曲下游 de novo 注释。）

与 gapfill 相同的招募设计（每个可填 gap 两翼 ~2 kb flank 建迷你
参考，minimap2 splice 预设比对），但证据是 cDNA：同时命中两翼的
转录本贡献其中间段。中间段接近零说明两翼是相邻外显子——gap
（基本）是内含子——默认 `--fill-mode cdna` 下直接对接（abut），
`--fill-mode exon-only` 下则保持开放。exon-only 还会跳过边缘
带有经典剪接位点（任一链上的 GT..AG）的 gap。

**参考位置引导的整基因放置。** 整个位点全是 N 的基因没有 flank 可供
转录本锚定，靠 flank 的填充永远无法恢复它——但转录组往往能把
这个基因完整组装出来。加 `--ref REF.fa` 后，rna-fill 对所有其他
机制都填不了的 gap 再做一遍放置：通过两翼 contig 的比对把 gap 括到
参考基因组上（bracket），落在括号内的参考基因与转录本（比对到参考）
匹配后，按参考顺序写进 gap——外显子为真实序列，内含子/基因间用
估计长度的 N 串占位（长度取自参考并截断封顶）。基因坐标来自
`--ref-gff`（若提供），否则用 miniprot 把 `--proteins` 贴回参考
基因组推导：

```bash
nearscaff rna-fill -a OUT/nearscaff.agp -q QUERY.fa \
    -T assembled.fa -o OUT_rnafill -t 16 \
    --internal-n 20 --fill-mode exon-only --one-sided \
    --ref REF.fa --proteins REF.pep.fa   # 或：--ref-gff REF.gff3
```

选项：`--fill-mode {cdna,exon-only}`；`--tx-preset`（默认
`splice:hq`，噪声大的 cDNA 用 `splice`）；`--tx-min-span INT`；
`--abut-window INT`（默认 30）；`--max-depth-factor FLOAT`（拒绝
多拷贝基因家族堆积转录本的 gap，默认 3.0，0 关闭）；
`--no-overlap-closure`；`--fill-tiers`（默认 `protein,asm5`）；
`--one-sided`（单侧延伸：转录本盖住一侧边缘并伸进 gap ≥500 bp 时
贴着该侧写出已知序列，残余 gap 保留）；`--internal-n LEN`（同时填充
组件内部 ≥LEN bp 的 N 串——短读组装的 scaffold 内部有大量估计长度的
N，HiFi 组装没有，默认 0 关闭）；`--ref REF.fa` + `--ref-gff GFF3` /
`--proteins`（参考位置引导的整基因放置；`--place-max-bracket`
默认 1 Mb，`--place-max-genes` 默认 50，`--place-min-ovlp`
默认 0.5；`--place-intact-cov` 默认 0.8——蛋白在 query 中已有
完整命中（无 X）的参考基因不重复放置，防旁系/挂载错位造重复；
`--place-max-spacer` 默认 50000，截断参考估计的内含子/基因间 N
占位长度）。有 `--ref` 时，长的 span/延伸填充还会对参考做校验：
≥500 bp 的填充若最优命中落在其 gap 的参考括号之外，判为疑似旁系
污染并驳回（该 gap 回落给放置机制；可用 `--no-fill-check` 关闭）。
需要 `PATH` 上有 `minimap2`。

输出：`nearscaff.rnafill.agp` 和 `nearscaff.rnafill.scaffolds.fa`
（填充 component 命名为 `*_gapfill_tx*`），另有按 span / abut /
overlap 闭合与内含子跳过分类的统计报告。

> **注意——填充物是 cDNA。** 默认 `cdna` 模式下填入的是真实的外
> 显子序列，但闭合 gap 内的内含子会被压缩丢掉。输出可用于提升注
> 释完整性（基因模型、BUSCO 蛋白、转录本回贴），但不要用于需要
> 真实内含子序列的分析。只动确定是外显子的 gap 时用
> `--fill-mode exon-only`。

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
- `intermediate/contig_alignments.tsv` —— Stage 0 与扩展轮缓存的逐
  contig 最佳比对，供最终精修复用、避免重复比对
- `intermediate/*.mmi` —— 可复用的、按预设区分的 minimap2 参考索引
- 使用 `--keep-intermediate` 时：`gene_anchors.tsv`、`synteny.blocks`

## 实际案例

两个 contig 级别的植物组装，使用同一套近缘染色体级别参考基因组
（14 条染色体）和参考蛋白集进行挂载，`--preset near`。

案例 1 —— 高度碎片化的组装：

| | 挂载前 | 挂载后 |
|---|---|---|
| 序列数 | 551,915 条 contig | 14 条染色体级 scaffold |
| 总大小 | 564.0 Mb | 挂载 469.8 Mb（83.3%） |
| N50 | 0.019 Mb | 35.5 Mb（scaffold） |
| BUSCO (embryophyta_odb10) | C:83.0% [S:80.7%, D:2.3%], F:13.3%, M:3.8% | C:95.1% [S:92.8%, D:2.4%], F:3.8%, M:1.1% |

案例 2 —— 长读长 contig 组装：

| | 挂载前 | 挂载后 |
|---|---|---|
| 序列数 | 17,594 条 contig | 14 条染色体级 scaffold |
| 总大小 | 456.0 Mb | 挂载 454.5 Mb（99.7%） |
| N50 | 2.95 Mb | 31.9 Mb（scaffold） |
| BUSCO (embryophyta_odb10) | C:97.7% [S:95.4%, D:2.3%], F:2.1%, M:0.2% | C:98.0% [S:95.6%, D:2.4%], F:1.9%, M:0.1% |

两个案例中 contig 都被挂载为与参考染色体一一对应的 scaffold，挂载后
完整 BUSCO 数量保持或有所提升。

## 与 RagTag 的对比

上述两个案例同时用 RagTag 基于同一参考基因组进行了挂载。BUSCO 评估
使用 embryophyta_odb10（n=1614，euk_genome_min 模式），只统计 14 条
染色体级 scaffold（未挂载的 contig 不计入）。

案例 1 —— 高度碎片化的组装：

| 指标 | RagTag | nearscaff |
|---|---|---|
| 挂载率（按碱基） | 58.4% | **83.3%** |
| scaffold N50 | 19 Mb | **35.5 Mb** |
| BUSCO — 染色体骨架 | C:86.6% [S:85.2%, D:1.4%], F:4.8%, M:8.6% | **C:95.1%** [S:92.8%, D:2.4%], F:3.8%, M:1.1% |

案例 2 —— 长读长 contig 组装：

| 指标 | RagTag | nearscaff |
|---|---|---|
| 挂载率（按碱基） | 87.2% | **99.7%** |
| scaffold N50 | 27 Mb | **31.9 Mb** |
| BUSCO — 染色体骨架 | C:97.2% [S:95.2%, D:2.0%], F:1.8%, M:1.0% | **C:98.0%** [S:95.6%, D:2.4%], F:1.9%, M:0.1% |

注：RagTag 有大量带基因的 contig 未被挂载（案例 1 中其染色体骨架
缺失 8.6% 的 BUSCO）；nearscaff 的核苷酸扩展轮次把这些 contig
大部分挂回了染色体。

### 挂载准确率

以独立的真值集评估挂载正确性：每个 query contig 取其对参考基因组的
最佳 minimap2 主比对（仅保留 mapq ≥ 30 且比对覆盖 ≥ 50% contig 长度
的记录）。指标为染色体一致率、方向一致率（在染色体一致的 contig 中）
和 scaffold 内顺序一致率（真值坐标的相邻对单调性）。

| 指标 | RagTag | nearscaff |
|---|---|---|
| 案例 1 —— 染色体一致率 | 99.86% | 98.83% |
| 案例 1 —— 方向一致率 | 99.92% | 93.73% |
| 案例 1 —— 顺序一致率 | 99.80% | 85.63% |
| 案例 2 —— 染色体一致率 | 97.35% | 89.96% |
| 案例 2 —— 方向一致率 | 99.74% | 96.62% |
| 案例 2 —— 顺序一致率 | 99.05% | 91.30% |

nearscaff 挂载的 contig 远多于 RagTag（案例 1：28.4 万 vs 4.5 万；
案例 2：1.7 万 vs 0.5 万）；多出的挂载来自 `asm5`–`asm20` 渐进扩展
轮次，用几个百分点的一致率换来了明显更高的挂载率和 BUSCO 完整度。
（上方一致率指标测自 0.3.0 版本的输出；挂载数量为当前版本数据。）

## 亚基因组分型（kmer-phase）

`nearscaff kmer-phase` 是独立于挂载流程的亚基因组分型模块，提供两种
工作模式（不随 `run` 默认运行）。

> **适用前提（请先读）**
>
> - **不推荐对碎片化组装做 contig 级分型。** 亚基因组的 k-mer 信号
>   只在 ≥ ~500 kb–1 Mb 的连续单元上才稳定存在；碎片化（contig
>   N50 仅几十 kb）的四倍体拿去分型，所有现有方法（包括
>   SubPhaser、Allo4D）都只能得到接近随机的结果，这是信号本身的
>   物理下限，不是参数问题。
> - 碎片化组装的正确路径是：**先用 `nearscaff run` 挂到染色体尺度，
>   再在染色体尺度上分型**（下面的模式 1）。
> - 模式 1（染色体尺度）要求组装已接近染色体级别；模式 2（碎片
>   级别，`--guide`）要求组装连续性尚可——至少要存在若干
>   ≥ 500 kb 的同源块（即 contig N50 在数百 kb 以上），且能提供
>   二倍体近缘参考。

**1. 染色体级别分型。** 适用于已挂到染色体级别的异源多倍体/混合组装，
最好配合同源染色体配对文件（每行一对同源染色体 ID）。核心算法是
同源组取向聚类：k-mer 的富集方向必须在**多数同源对**中保持一致才被
采纳，谱方法估计全局取向，bootstrap 给出逐条置信度，同源对被强制
拆分、弱置信成员如实拒判。

```bash
nearscaff kmer-phase -q assembly.fa -c homology.cfg -o phasing.tsv -t 16
```

**2. 碎片级别分型（块引导）。** 适用于 contig 级别的多倍体/混合组装，
需要一个二倍体近缘参考（基因组 + 蛋白）。全自动链条：miniprot 注释
→ jcvi 共线同源块 → 块级取向聚类（给长 contig 打上近乎完美的标签）
→ EM LLR 精修，把信号传播到没有自有信号的短 contig 上。

```bash
nearscaff kmer-phase -q mixed.fa \
    --guide diploid.pep --guide-ref diploid.genome.fa \
    -o phasing.tsv --outdir guide_work -t 16
```

依赖：jellyfish、numpy/scipy/sklearn（聚类）；块引导模式另需
miniprot、gffread 和 jcvi（含 LAST，`NEARSCAFF_JCVI` 可指定命令前缀，
其他工具可用 `NEARSCAFF_<工具名>` 覆盖）。

### 分型准确率对比

三个匿名评测集（伪异源四倍体，真值 = contig/染色体来源的物种）：
案例 A 为近缘组合（碎片化）、案例 B 为中等分化组合（碎片化）、
案例 C 为三代长读长组合。

染色体尺度（correct / total）：

| 方法 | 案例 A（12 条） | 案例 B（16 条） |
|---|---|---|
| SubPhaser | 7/12（11:1 塌缩） | 15/16 |
| kmer-phase | **12/12** | **15/15**（1 条低置信拒判） |

碎片（contig）尺度（真值 = 每条 contig 的来源物种）：

| 方法 | 案例 A | 案例 B | 案例 C |
|---|---|---|---|
| SubPhaser | 完全塌缩，无法分型 | 完全塌缩，无法分型 | — |
| Allo4D | 正确率 55%，且只覆盖 ~30% 的基因组 | 正确率 54%，只覆盖 17% | 正确率 74%，只覆盖 48% |
| kmer-phase --guide | **68% 的 contig 分对；按碱基算 86% 分对** | 不适用¹ | **98% 的 contig 分对；按碱基算 99.95% 分对** |

「按碱基算」指长度加权准确率——越长的 contig 权重越大。它比
contig 准确率高，说明分错的主要是短 contig：案例 A 中基因组绝大
部分碱基（86%）已归属正确，剩余错误集中在可信度本就该存疑的
碎片上；案例 C 则基本全对。

¹ 碎片化混样中不存在 ≥500 kb 的可分辨同源块，信号在物理上不足；
建议先用 nearscaff 挂到染色体尺度再分型（案例 B 染色体尺度 15/15）。

碎片尺度分型的效果由组装连续性决定：≥500 kb 的可分辨同源块越多、
覆盖越广，全基因组准确率越高（块 ≥500 kb 时块级准确率即达 97%）。
置信度低的输出会被如实标注或拒判。

## 更新日志

### 未发布

- **`rna-fill`：转录组引导的基因区 gap 填补。** 由 gap 两翼招募
  转录本（短读组装或现成转录本）；支持 span 闭合、剪接感知的
  abut/overlap 处理、`--one-sided` 单侧延伸、`--internal-n`（炸开
  contig 内部 N 串）和 `--fill-mode {cdna,exon-only}`。
- **参考位置引导的整基因放置（`rna-fill --ref`）。** flank 机制
  无法触达的 gap（整个基因座全是 N）借两翼 contig 的比对括到参考
  基因组上；括号内的参考基因（`--ref-gff`，或用 `--proteins` 经
  miniprot 推导）匹配转录本后按参考顺序写入——真实外显子序列 +
  带 GT..AG 剪接边的内含子占位。intact 反重复检查
  （`--place-intact-cov`）跳过 query 中已存在的基因；其余用
  `--place-max-bracket/-genes/-min-ovlp/-max-spacer` 调节。
- **fill-check 后置校验（`--ref` 时默认开启）。** ≥500 bp 的
  span/延伸填充会对参考做校验：最优命中落在 gap 参考括号之外
  （±50 kb）的判为疑似旁系污染并驳回，gap 回落给放置机制
  （`--no-fill-check` 关闭）。
- **`gapfill`：长读 gap 闭合**（`--method sr/lr/endjoin`）。
- **`--keep-unplaced [MINLEN]`。** 从未进图的 query contig 以原名
  singleton 追加进 AGP（默认关闭）——高度碎片化的 query 否则会
  从输出中丢失相当比例的组装序列。

### 0.3.1

Stage 1（核苷酸扩展）的 CPU/资源优化——在挂载精度不降（实测提升）
的前提下大幅缩短运行时间：

- 每轮扩展改为对可复用的、按预设区分的 minimap2 索引做一次全参考
  比对，取代旧的「逐区域 × 逐轮」比对循环（旧路径在碎片化组装上
  是病态的）。
- Stage 0 的 contig→参考比对结果写入缓存
  （`intermediate/contig_alignments.tsv`），最终精修直接复用；只对
  无缓存的 contig 重新比对（并用 `--secondary=no`）。
- 扩展轮保留至多 `-N 5` 个次级比对（`--secondary-alignments`），
  且添加扩展边时考虑全部比对——主比对离所有 scaffold 都远的 contig
  也能通过次级比对挂回。
- `--nucleotide-passes` 默认改为 `asm5 asm20`（`asm10` 轮冗余）。
- 修复一次性用 `samtools faidx` 抽取极多 contig 时参数列表过长导致
  的崩溃（现在回退为流式扫描）。

在下面两个实际案例中实测（同一台机器）：Stage 1 从「跑不完 /
数小时」降到约 32 秒（案例 2）和约 574 秒（案例 1），染色体尺度
结构一致（14 条 scaffold），BUSCO 完整度持平或更好：

| 案例 | BUSCO — 0.3.0 | BUSCO — 0.3.1 |
|---|---|---|
| 案例 1（碎片化组装） | C:92.1% [S:90.1%, D:2.0%], F:5.4%, M:2.5% | **C:95.1%** [S:92.8%, D:2.4%], F:3.8%, M:1.1% |
| 案例 2（长读长组装） | C:97.6% [S:95.5%, D:2.0%], F:1.7%, M:0.7% | **C:98.0%** [S:95.6%, D:2.4%], F:1.9%, M:0.1% |

## 开发

运行测试（当 miniprot/minimap2/samtools 不在 PATH 中时，集成测试会
自动跳过）：

```bash
python -m pytest tests -q
```
