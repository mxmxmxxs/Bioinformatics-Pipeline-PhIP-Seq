
#PhIP pipeline: Build bowtie2 index, run alignment, convert SAM to BAM, and generate idxstats

threads=8 #number of threads to use, can be adjusted

work_dir="/Users/mima/Desktop/Analysis_scripts_and_input_files/" #working directory for files
mkdir -p $work_dir #make working directory if it doesn't exist
bw2_base=$work_dir"vir3" #bowtie2 index base path, include / if no seperator in work_dir
fasta="/Users/mima/Desktop/Analysis_scripts_and_input_files/vir3.fasta" #path to fasta file for building bowtie2 index

bowtie2-build --threads $threads $fasta $bw2_base #build bowtie2 index from fasta file

echo "Bowtie2 index base path:"$bw2_base #print bowtie2 index base path to see where in the process we are while running the script

fastq_dir="/Users/mima/Desktop/VirScan/VirScan_pipeline/Example_data/fastq_files/" #directory containing fastq files

output_dir="/Users/mima/Desktop/Masterthesis/000_SpecialCourse_KU/Code/data/" #directory to save output files

for fastq_file in "$fastq_dir"*_R1.fastq; do #loop through all fastq files in the fastq directory
    filename=$(basename "$fastq_file") 
    sample_name=${filename%%_R1.fastq} #extract sample name from filename
    echo "Sample: $sample_name" #print sample name to see progress

    file_path=$fastq_dir/$sample_name #base file path without extension

    alignment=$work_dir"alignment.sam" #path to alignment file
    echo "Alignment file:"$alignment #print alignment file path to see progress
    report=$file_path".txt" #report file path

    samfile=$alignment #sam file path
    bamfile=$file_path"_unsorted.bam" #unsorted bam file path
    sorted_file=$file_path"_sorted.bam" #sorted bam file path
    
    bowtie2 --trim3 25 --very-sensitive -p $threads -x $bw2_base -U $fastq_file -S $alignment        #trim N bases from 3' end of each read before alignment, run bowtie2 alignment
    
    echo "Bowtie2 mapping completed for $sample_name" #print progress message for this sample

    samtools view -bS -@ $threads $samfile > $bamfile #convert sam to bam
    samtools sort -o $sorted_file -@ $threads $bamfile #sort bam file
    samtools index -b $sorted_file #index sorted bam file
    samtools idxstats -@ $threads $sorted_file | sed 's/\t/;/g' > $output_dir$sample_name"_idxstats.txt" #generate idxstats and save to output directory

    echo "SAM to BAM conversion, sorting, and indexing completed for $sample_name" #print progress message for this sample
    
done

echo "all done, time for analysis" #print final message when all samples are done

#Run analysis and visualization
echo "Starting analysis and visualization" #print starting message
python3 /Users/mima/Desktop/Masterthesis/000_SpecialCourse_KU/Code/analysis_and_Visualization.py

echo "Analysis and visualization complete!" #print final completion message
