''' Contains numerous common tools for bioinformatics analysis


Functions within:

rev_comp -- return reverse complement of DNA sequence

fasta_fixer -- create a reformated fasta files which contain return characters within the sequences

fasta_from_IDs -- filter a fasta file to contain only sequences within the input ID list

alignment_gap_makser -- mask gapped positions in a MSA

sam_to_fasta -- convert a .sam file to .fasta

'''




from sys import argv





aa_list = ['A','R','N','D','C','E','Q','G','H','I','L','K','M','F','P','S','T','W','Y','V']









def rev_comp(input_seq):
	rev_comp_dict = {'A':'T','C':'G','G':'C','T':'A'}
	return ''.join(rev_comp_dict[nuc] for nuc in input_seq[::-1])
	










def fasta_fixer(input_file,output_file):
	'''Reformat fasta file to remove newline characters within sequence

	Keyword arguments:
		input -- fasta file to fix location
		output -- fixed fasta file location
	
	'''
	out = open(output_file,'w')

	for i,l in enumerate(open(input_file,'U')):
		if l[0] == '>':
			if i == 0:
				out.write(l)
			else:
				out.write('\n'+l)
		else:
			out.write(l.strip())
	out.close()










def fasta_from_IDs(ID_file,fasta_file,output_file):
	
	'''Filter fasta file for sequences within input ID list
	
	'''
	
	ID_list = []
	for l in open(ID_file,'U'):
		ID_list.append(l.rstrip())
	
	out = open(output_file,'w')
	
	keep = False
	for l in open(fasta_file,'U'):
		if l[0] == '>':
			header = l
			name = l.replace('>','')
			name = name.rstrip()
			if name in ID_list:
				keep = True
		else:
			if keep == True:
				out.write(str(header) + str(l))
				print('yes')
			keep = False
	out.close()











def alignment_gap_masker (input_file,output_file,threshold):
    out = open(output_file,'w')

    input_lines = 0
    input_alignment_length = 0
    for l in open(input_file,'U'):
        if input_lines == 1:
            input_alignment_length = len(l.strip())
        input_lines += 1
    keep_string = ''
    seq_strings = ['' for x in range(input_alignment_length)]
    seq_names = []

    for counter_1,l in enumerate(open(input_file,'U')):
        if l[0] == '>':
            seq_names.append(l)
        else:
            l_stripped = l.strip()
            for i,letter in enumerate(l_stripped):
            	seq_strings[i] += letter
        if counter_1 % 1000 == 0:
            print('Step 1 is ' + str((counter_1/float(input_lines))*100) + '% complete')

    print('Step 1 done')            

    for counter_2,column in enumerate(seq_strings):
        if column.count('-')/float(len(column)) > float(threshold):
            keep_string += 'F'
        else:
            keep_string += 'T'
        if counter_2 % 500 == 0:
            print('Step 2 is  ' + str((counter_2/float(input_alignment_length))*100) + '% complete')

    print('Step 2 done')
        
    for counter_3,seq in enumerate(open(input_file,'U')):
        if seq[0] == '>':
            out.write(seq)
        else:
            seq_stripped = seq.strip()
            for pos,aa in enumerate(seq_stripped):
                if keep_string[pos] == 'T':
                    out.write(aa)
            out.write('\n')
        if counter_3 % 1000 == 0:
            print('Step 3 is ' + str((counter_3/float(input_lines))*100) + '% complete')
    out.close()











def sam_to_fasta(input_sam,output_fasta):
	out = open(output_fasta,'w')
	for l in open(input_sam,'U'):
		if l.lstrip()[0] != '@':
			l = l.split('\t')
			out.write('>'+l[0]+'\n')
			out.write(l[9]+'\n')
	out.close()









if __name__ == '__main__':
	import subprocess
	fasta_from_IDs('./bioinfo_tools_tests/test_IDs.txt', './bioinfo_tools_tests/test_fasta.fasta', './bioinfo_tools_tests/test_output.fasta')
	print('\nfasta_from_ID test results:\n\n')
	subprocess.call('head ./bioinfo_tools_tests/test_output.fasta', shell=True)
	print('\n'+'-'*15)