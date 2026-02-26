### 介绍
该工程用于自连序列的拆分，自连序列中的每段序列为primer序列。输入相应待拆分文件和primer文件。
比如：
ATC-------CGGA 或者是 TCCG-------GAT
首尾两端都需要包含primer序列
primer序列保存在一个单独的fasta文件中，分别以"_L" 和 "_T"结尾，
其中_L是上游的序列，_T为下游的序列，注意与barcode当中不同的是，_T为
反义互补链的5'->3'的序列，在barcode文件中_R为正链的5'->3'的序列。
>P01_bnlg439w1_L <br>
AGTTGACATCGCCATCTTGGTGAC <br>
>P01_bnlg439w1_T <br>
GAACAAGCCCTTAGCGGGTTGTC <br>

### 参数
>
--pipeline_version  default_value ="1"， 内部拆分的核心逻辑版本 \
--primer 指定 primer 序列文件 (FASTA 格式) \
--max_distance 匹配pattern时允许的最大错误率 \
--input_file 指定输入序列文件 (FASTA/FASTQ/BAM) \
--output_folder 输出文件目录 \
--log_folder log文件的输出目录 \
--keep_primer 是否将primer序列保存在输出序列中，通常都需要开启 \
--tail_cutoff 由于primer自连序列可能在自连处存在连接问题，可以剔除掉若干碱基，来增加找到的序列的数量。设置的时候
要注意输入序列的质量，以及primer之间的相似度问题。 \
--reservesed_threads 保留的线程数量 \
--output_format 输出的数据格式，支持 "fasta", "fastq", "bam" \
--min_subread_len 最小子读长度，默认50,如果是保留primer，是包含在里面的，需要适量增大该参数
>