#!/usr/bin/env python2

## ----------------------IMPORT NECESSARY PACKAGES----------------------- 
from optparse import OptionParser
import os
from os import popen
import shutil
import subprocess
import multiprocessing
import time
start_time=time.time()

def Translate(input_dir,output_dir):
	make_temp="mkdir "+output_dir+"/translate"
	os.system(make_temp)
	Translate_output_path=output_dir+"/translate"
	for filename in os.listdir(input_dir):
		sample=filename.split("_")[0]
		filepath=input_dir+"/"+filename
    		file = open(filepath,mode='r')
    		all_of_it=file.read()
    		file.close
    		os.system('curl -s -d "dna_sequence='+str(all_of_it)+'&output_format=fasta" https://web.expasy.org/cgi-bin/translate/dna2aa.cgi > '+Translate_output_path+'/'+sample+'_proteinconv.fasta')

##---------------------RUN CARD------------------------------
def CARD(input_dir,output_dir):
	print("Running CARD-rgi......")
	make_temp="mkdir "+output_dir+"/CARD"
	os.system(make_temp)
	CARD_output_path=output_dir+"/CARD"
	hold_files=os.listdir(input_dir)
	multi_list=[hold_files[x:x+4] for x in range(0,len(hold_files),4)]
	for x in multi_list:
		processes=[]
		for files in x:
			prefix=files.split("_")[0]
			input_path=input_dir+"/"+files
			final_output=CARD_output_path+"/"+prefix+"_CARD-RGI_coding"
			command = ["rgi","-i",input_path,"-o",final_output]
			f=os.tmpfile()
			p=subprocess.Popen(command,stdout=f)
			processes.append((p,f))
		for p,f in processes:
			p.wait()
			f.seek(0)
			f.close()
	for file in os.listdir(CARD_output_path):
		file_path=CARD_output_path+'/'+file
		if os.stat(file_path).st_size == 0:
			command="rm "+file_path
			os.system(command)
	return(CARD_output_path)

##-----------------------RUN VFDB--------------------------------
def VFDB(input_dir,output_dir):
	print("Running blastp on VFDB......")
	make_temp="mkdir "+output_dir+"/VFDB"
	os.system(make_temp)
	VFDB_output_path=output_dir+"/VFDB"
	hold_files=os.listdir(input_dir)
	multi_list=[hold_files[x:x+4] for x in range(0,len(hold_files),4)]
	
	for x in multi_list:
		processes=[]
		for files in x:
			prefix=files.split("_")[0]
			command = ["blastn", "-db", "../Tools/VFDB/Virulence_Factors_core", "-query",input_dir+"/"+files,"-out",VFDB_output_path+"/"+prefix+"_VFDB_coding", "-max_hsps","1","-max_target_seqs","1","-num_threads","4","-evalue","1e-5"]
			f=os.tmpfile()
			p=subprocess.Popen(command,stdout=f)
			processes.append((p,f))
		for p,f  in processes:
			p.wait()
			f.seek(0)
			f.close()
	for file in os.listdir(VFDB_output_path):
		file_path=VFDB_output_path+"/"+file
		if os.stat(file_path).st_size == 0:
			command="rm "+file_path
			os.system(command)
	return(VFDB_output_path)

##--------------------RUN PILERCR----------------------------------
def Pilercr(input_dir,output_dir):
	print("Running Piler-CR.....")
	make_temp="mkdir "+output_dir+"/pilercr"
	os.system(make_temp)
	pilercr_output_path=output_dir+"/pilercr"
	hold_files=os.listdir(input_dir)
	multi_list=[hold_files[x:x+4] for x in range(0,len(hold_files),4)]
	
	for x in multi_list:
		processes=[]
		for filename in x:
			prefix = filename.split("_")[0]
			file_path=input_dir+"/"+filename
			command = ["../Tools/pilercr/pilercr1.06/pilercr", "-in",file_path, "-out",pilercr_output_path+"/"+prefix+"_pilercr_coding"]
			f=os.tmpfile()
			p=subprocess.Popen(command,stdout=f)
			processes.append((p,f))
		for p,f in processes:
			p.wait()
			f.seek(0)
			f.close()
	for file in os.listdir(pilercr_output_path):
		file_path=pilercr_output_path+"/"+file
		if os.stat(file_path).st_size == 0:
			command="rm "+file_path
			os.system(command)
	return(pilercr_output_path)

##-------------------------FORMAT OUTPUT FROM UCLUST SO IT CAN BE EASILY SPLIT FOR MULTIPROCESSING WITH EGGNOG----------------------------
def format_uclust(input,output_dir):
	file_split=open(input,"r")
	new_file=open(output_dir+"/USEARCH/All_format.txt","w+")
	file_split_read=file_split.readlines()
	count=0
	for lines in file_split_read:
		if not ">" in lines:
			lines_reformat=lines.rstrip()
			new_file.write(lines_reformat)
		elif ">" in lines:
			if count ==0:
				new_file.write(lines)
			else:
				output="\n"+lines
				new_file.write(output)
			count=1
	new_file.close()
	file_split.close()
	return(output_dir+"/USEARCH/All_format.txt")

##----------------------RUN EGGNOG----------------------------------------------
def eggnog(input,output_dir):
	splitLen=100
	at=1
	make_temp="mkdir "+output_dir+"/eggNOG"
        os.system(make_temp)
	make_temp="mkdir "+output_dir+"/eggNOG_split"
	os.system(make_temp)
	outputBase=output_dir+"/eggNOG_split/eggNOG_split"
	eggNOG_annotation_path=output_dir+"/eggNOG"
	input_split=open(input,"r").read().split("\n")
	for lines in range(0,len(input_split),splitLen):
		outputData=input_split[lines:lines+splitLen]
		output=open(outputBase+str(at)+'.txt','w')
		output.write('\n'.join(outputData))
		output.close()
		at +=1
	hold_files=os.listdir(output_dir+"/eggNOG_split")
	multi_list=[hold_files[x:x+4] for x in range(0,len(hold_files),4)]
	for x in multi_list:
		processes=[]
		for files in x:
			prefix=files.split("_")[1]
			file_path=output_dir+"/eggNOG_split/"+files
			call_eggNog=["emapper.py", "-i",file_path,"--data_dir", "../Tools/eggnog_db","-m","diamond","--dmnd_db","/home/projects/group-c/Team3-FunctionalAnnotation/Tools/eggnog_db/eggnog_proteins.dmnd","--translate", "-d", "bact", "--output_dir",eggNOG_annotation_path, "-o","eggNOG_"+prefix]
			f=os.tmpfile()
			p=subprocess.Popen(call_eggNog,stdout=f)
			processes.append((p,f))
		for p,f in processes:
			p.wait()
			f.seek(0)
			f.close()
	
	command="for i in "+output_dir+"/eggNOG_split; do cat $i > "+output_dir+"/eggNOG/eggNOG_combined.txt ; done"
	os.system(commnand)
	return(output_dir+'/eggNOG/eggNOG_combined.txt')

##---------------------------RUN SIGNALP----------------------------------
def SignalP(input_dir,output_dir):
	make_temp="mkdir "+output_dir+"/SignalP"
	os.system(make_temp)
	SignalP_path=output_dir+"/signalP"
	for filename in os.listdir(input_dir):
		sample=filename.split("_")[0]
    		os.system("/home/projects/group-c/Team3-FunctionalAnnotation/Tools/signalp-5.0b/bin/signalp -fasta "+input_dir+"/"+filename+" -org gram+ -format short -gff3 -prefix "+output_dir+"/SignalP/"+sample+"_signalp")

##--------------------------RUN UCLUST FOR EGGNOG INPUT---------------------------
def uclust(input,output_dir):
	make_temp="mkdir "+output_dir+"/USEARCH"
	os.system(make_temp)

	combine_path=output_dir+"/USEARCH/All.fasta"
	output=open(combine_path,"w+")
	for filename in os.listdir(input):
		filepath=input+"/"+filename
		genes_file=open(filepath,"r")
		sample=filename.split("_")
		sample=sample[0]
        	for line in genes_file:
			if ">" in line:
				labled_line=line.rstrip()
				labled_line=labled_line+"_"+sample+"\n"
                        	output.write(labled_line)
                	else:
                     		output.write(line)
        	genes_file.close()
	output.close()
	cluster_path=output_dir+"/USEARCH/All_clustered.fasta"
	uc_path=output_dir+"/USEARCH/clusters.uc"
	call_usearch="../Tools/USEARCH/usearch -cluster_fast "+combine_path+" -id 0.97 -centroids "+cluster_path+" -uc "+uc_path 
	os.system(call_usearch)
	return(cluster_path,uc_path)

##--------------MAP NODES OF CLUSTERS FOR LATER MAPPING OF EGGNOG OUTPUTS TO EACH FILE---------------------
def mapNodes(input_clusters):
	path_to_clusters=input_clusters
	coding_cluster=open(path_to_clusters,"r")
	coding_cluster_read=coding_cluster.readlines()
	cluster_dict={}
	for line_cluster in coding_cluster_read:
		hold_cluster=line_cluster.split("\t")
		node_cluster=hold_cluster[9]
		if '*' in node_cluster:
			node_cluster=hold_cluster[8]
		node_cluster=node_cluster.rstrip()
		node_cluster=node_cluster.split("_")
		node_cluster=node_cluster[:-1]
		join_node="_"
		node_cluster=join_node.join(node_cluster)
		if node_cluster in cluster_dict:
			cluster_dict[node_cluster]=(cluster_dict.get(node_cluster))+([line_cluster])
		else:
			cluster_dict[node_cluster]=[line_cluster]
	coding_cluster.close()
	return(cluster_dict)

##---------------FORMAT OUTPUT OF PILERCR SO THAT IT IS IN GFF FORMAT---------------------------------
def formatPilercr(input_dir,output_dir):
	print("formatting Pilercr....")
	make_temp="mkdir -p "+output_dir+"/format/Pilercr"
	os.system(make_temp)
	path_to_Pilercr=input_dir
	path_to_output=output_dir+"/format/Pilercr"
	for filename in os.listdir(path_to_Pilercr):
		sample=filename.split("_")[0]
		path_to_file=path_to_Pilercr+"/"+filename
		Pilercr_file=open(path_to_file,"r")
		Pilercr_file_read=Pilercr_file.readlines()
		Pilercr_output_path=path_to_output+"/"+sample+"_pilercr.gff"
		Pilercr_output=open(Pilercr_output_path,"w+")
		hold_nodes=[]
		for line in range(len(Pilercr_file_read)):
			if Pilercr_file_read[line].startswith("Array ") and ">" in Pilercr_file_read[line+1]:
				node=Pilercr_file_read[line+1]
				node_full=node.split("\t")[0]
				hold_nodes=hold_nodes+[node_full]
			if 'SUMMARY BY POSITION' in Pilercr_file_read[line]:
				index=0
				while line+index < len(Pilercr_file_read):
					for nodes in hold_nodes:
					
						if nodes in Pilercr_file_read[line+index]:
							node=nodes
							info_line=Pilercr_file_read[line+index+4].split(" ")
							hold_digit=[]
							for i in info_line:
								if i.isdigit():
									hold_digit=hold_digit+[i]
							start=hold_digit[1]
							end=str(int(hold_digit[2])+int(start))
							gff_line=node.rstrip()+"\tPilercr\tCRISPR_array\t"+start+"\t"+end+"\t.\t.\t.\tfeature=Putative CRISPR array"
							print(gff_line)
							Pilercr_output.write(gff_line)
					index=index+1
		Pilercr_output.close()
		Pilercr_file.close()

##-------------------------FORMAT OUTPUT OF EGGNOG SO IT IS IN GFF FORMAT--------------------------
def formateggNOG(input,cluster_dict,output_dir):
	make_temp="mkdir -p "+output_dir+"/format/eggNOG"
	eggNog_file=open(input,"r")
	os.system(make_temp)
	path_to_eggnog_out=input
	for line in eggNog_file:
		if line.startswith("NODE"):
			annotation_split=line.split("\t")
			node_eggnog=annotation_split[0]			
			node_eggnog=node_eggnog.split("_")
			node_eggnog=node_eggnog[:-1]
			if len(node_eggnog) > 5:
				join_node=("_")
				node_eggnog=join_node.join(node_eggnog)
				if node_eggnog in cluster_dict:
					loop_lines=cluster_dict.get(node_eggnog)
					for line2 in loop_lines:
						cluster_split=line2.split("\t")
						if node_eggnog in cluster_split[9] or node_eggnog in cluster_split[8]:
							mapped_node=cluster_split[8]
							mapped_node=mapped_node.split("_")
							sample_node=mapped_node[-1]
							mapped_node=mapped_node[:-1]
							join_node=("_")
							mapped_node=join_node.join(mapped_node)
							gene_start_end=mapped_node.split(":")
							gene_start_end=(gene_start_end[1]).split("-")
							start="1"
							end=str(int(gene_start_end[1])-int(gene_start_end[0]))
							annotation=annotation_split[1:]
							join_annotation=";"
							annotation=join_annotation.join(annotation)
							output_path=output_dir+"/format/eggNOG/"+sample_node+"_eggNOG.gff"
							if os.path.exists(output_path):
								gff_eggnog=open(output_path,"a")
							else:
								gff_eggnog=open(output_path,"w+")
							gff_write=mapped_node+"\teggNOG\tprotein_match\t"+start+"\t"+end+"\t.\t.\t.\t"+annotation
							#print(gff_write)
							gff_eggnog.write(gff_write)
							gff_eggnog.close()
	command="for i in "+output_dir+"/format/eggNOG/*; do sort -u $i > ${i}_sorted;done"
	command2="rm "+output_dir+"/format/eggNOG/*.gff"
	os.system(command)
	os.system(command2)

##-------------------FORMAT CARD OUTPUT SO IT IS IN ONE LINE GFF FORMAT-------------------------------
def formatCARD(input_dir,output_dir):
	make_temp="mkdir -p "+output_dir+"/format/CARD"
	os.system(make_temp)
	path_to_CARD=input_dir
	print("formatting...card")
	for filename in os.listdir(path_to_CARD):
		if "gff" in filename:
			sample=filename.split("_")[0]
			path_CARD_file=path_to_CARD+"/"+filename
			CARD_file=open(path_CARD_file,"r")
			CARD_output_path_gff=output_dir+"/format/CARD/"+sample+"_card_mapped.gff"
			CARD_output_gff=open(CARD_output_path_gff,"w+")
			for line in CARD_file:
				if line.startswith("NODE"):
					annotation=line.split("\t")
					node=annotation[0]
					node=node.split("_")
					node=node[:-1]
					node_join="_"
					node=node_join.join(node)
					annotation[0]=node
					annotation_join="\t"
					#print(annotation)
					annotation=annotation_join.join(annotation)
					gff_line=annotation
					#print(gff_line)
					CARD_output_gff.write(gff_line)
			CARD_output_gff.close()	

##-------------------FORMAT VFDB SO IT IS IN GFF FORMAT-------------------------------------	
def formatVFDB(input_dir,output_dir):
	VFDB_dir=input_dir
	make_temp="mkdir -p "+output_dir+"/format/VFDB"
	os.system(make_temp)
	print("formatting...VFDB")
	for filename in os.listdir(VFDB_dir):
		full_file_path=VFDB_dir+"/"+filename
		VFDB_file=open(full_file_path,"r")
		VFDB_file_read=VFDB_file.readlines()
		sample=filename.split("_")[0]
		output_path_gff=output_dir+"/format/VFDB/"+sample+"_VFDB.gff"
		gff_VFDB=open(output_path_gff,"w+") 
		for index in range(len(VFDB_file_read)):
			if "Query=" in VFDB_file_read[index] and "Sequences producing" in VFDB_file_read[index+4]:
				node=VFDB_file_read[index].split(" ")
				node=node[1]
				count=9
				start_bool=False
				while not "Lambda" in VFDB_file_read[index+count]:
					if "Score =" in VFDB_file_read[index+count]:
						score=VFDB_file_read[index+count]
						score=score.split(" ")
						score=score[-1].rstrip()
						strand=VFDB_file_read[index+count+2]
						strand=strand.split("/")[1].rstrip()
						if strand=="Plus":
							strand="+"
						if strand=="Minus":
							strand="-"
					if VFDB_file_read[index+count].startswith("Query"):
						split_line=VFDB_file_read[index+count].split(" ")
						if not start_bool:
							start_bool=True
							start=split_line[2]
						end=split_line[-1]
						end=end.rstrip()
					count=count+1
				gff_write=node.rstrip()+"\tVFDB\tprotein_match\t"+start+"\t"+end+"\t"+score+"\t"+strand+"\t.\tstitle="+VFDB_file_read[index+9].lstrip('>')
				gff_VFDB.write(gff_write)
		
		VFDB_file.close()
		gff_VFDB.close()

##-----------------------------------FORMAT SIGNALP SO IT IS IN GFF FORMAT------------------------------
def formatSignalP(input_dir,output_dir):
	print("Formatting...SignalP")
	path_to_SignalP=input_dir
	make_temp="mkdir -p "+output_dir+"/format/SignalP"
	os.system(make_temp)
	for filename in os.listdir(path_to_SignalP):
		if "gff" in filename:
			path_SignalP_file=path_to_SignalP+"/"+filename
			SignalP_file=open(path_SignalP_file,"r")
			sample_name=filename.split("_")
			sample_name=sample_name[1]
			SignalP_gff_path=output_dir+"/format/SignalP/"+sample_name+"_SignalP.gff"
			SignalP_gff=open(SignalP_gff_path,"w+")
			for line in SignalP_file:
				if line.startswith("NODE"):
					SignalP_gff.write(line)
			SignalP_gff.close()
			SignalP_file.close()

##-------------------------------MERGE ALL GFF FORMATTED FILES INTO ONE SO THAT NODES ARE IN ORDER FOR EACH SAMPLE------------------------
def mergeGff(output_dir):
	print("Merging Files.....")
	make_temp="mkdir -p ./Outputs/merged"
	os.system(make_temp)
	sample_dict={}
	mapped_dir_path=output_dir+"/format"
	for dir_name in os.listdir(mapped_dir_path):
		if dir_name=="eggNog" or dir_name=="CARD" or dir_name=="VFDB" or dir_name=="eggNOG" or dir_name=="SignalP" or dir_name=="Pilercr":
			tool_dir_path=mapped_dir_path+"/"+dir_name
			for filename in os.listdir(tool_dir_path):
				if "gff" in filename and not "getorf" in filename:
					filename_split=filename.split("_")
					sample=filename_split[0]
					tool=filename_split[1]
					file_path=tool_dir_path+"/"+filename
					annotation_file=open(file_path,"r")
					annotation_file_read=annotation_file.readlines()
					if sample in sample_dict:
						feature_dict=sample_dict.get(sample)
						for line in annotation_file_read:
							line_split=line.split("\t")
							node=line_split[0]
							write_dict=line
							if node in feature_dict:
								feature_dict[node]=feature_dict.get(node)+[write_dict]
							else:
								feature_dict[node]=[write_dict]
					if sample not in sample_dict:
						feature_dict={}
						for line in annotation_file_read:
							line_split=line.split("\t")
							node=line_split[0]
							#print(node)
							write_dict=line
							if node in feature_dict:
								feature_dict[node]=feature_dict.get(node)+[write_dict]
							else:
								feature_dict[node]=[write_dict]
						sample_dict[sample]=feature_dict
					annotation_file.close()


	for key,value in sample_dict.items():
		sorted_dict={}
		feature_dict=value
		output_path=("./Outputs/merged/"+key+"_merged.gff")
		output_file=open(output_path,"w+")
		sorted_list=sorted(feature_dict, key=lambda key1: int(key1.split("_")[1]))
		for i in sorted_list:
			node=i.split("_")[0:2]
			join_node="_"
			node=join_node.join(node)		
			if node in sorted_dict:
				sorted_dict[node].append(i)
			else:
				sorted_dict[node]=[i]
		sorted_list=sorted(sorted_dict, key=lambda key1: int(key1.split("_")[1]))
		for value_item in sorted_list:
			item=sorted_dict.get(value_item)
			#print(item)
			value_sort=sorted(item, key=lambda key1: (int((key1.split(":")[1]).split("-")[1])))
			for get_item in value_sort:
				loop_lines=feature_dict.get(get_item)
				node=get_item
				header="##sequence_region "+node+"\n"
				output_file.write(header)
				for annotation in loop_lines:
					output_file.write(annotation)
		output_file.close()
		make_temp="mkdir -p ./Outputs/CARD"
		os.system(make_temp)
		make_temp="mkdir -p ./Outputs/VFDB"
		os.system(make_temp)
		copy_card="for i in "+output_dir+"/format/CARD/*; do cp $i ./Outputs/CARD; done"
		#print(copy_card)
		os.system(copy_card)
		copy_vfdb="for i in "+output_dir+"/format/VFDB/*; do cp $i ./Outputs/VFDB; done"
		os.system(copy_vfdb)
		#print(copy_vfdb)
def GetFASTA(input,dir,toolname):
	make_temp="mkdir -p ./Outputs/fasta/"+toolname
	os.system(make_temp)
	for filename in os.listdir(dir):
		if "gff" in filename:
			gff_file_path=dir+"/"+filename
			print(filename)
			sample=filename.split("_")[0]
			for file in os.listdir(input):
				if sample in file:
					fasta_file_path=input+"/"+file
					output_path="./Outputs/fasta/"+toolname+"/"+sample+"_"+toolname+".fna"
					command="bedtools getfasta -fi "+fasta_file_path+" -bed "+gff_file_path+" > "+output_path
			#print(command)
					os.system(command)
	command="rm "+input+"/*.fai"
	os.system(command)
def GetFASTA_all(input,dir):
	make_temp="mkdir -p ./Outputs/fasta/merged"
	os.system(make_temp)

	for filename in os.listdir(dir):
		sample=filename.split("_")[0]
		path_to_file=dir+"/"+filename
		for file in os.listdir(input):
			if sample in file:
				path_to_fna=input+"/"+file
				fna_file=open(path_to_fna,"r")
				fna_file_read=fna_file.readlines()
				output_path="./Outputs/fasta/merged/"+sample+"_for_comparative.fna"
				output_file=open(output_path,"w+")
				gff_file=open(path_to_file,"r")
				for line in gff_file:
					if "#" in line:
						node=line.split(" ")[1]
						for line2 in range(len(fna_file_read)):
							if node in fna_file_read[line2]:
								node_line=str(fna_file_read[line2])
								sequence_line=str(fna_file_read[line2+1])
								output_file.write(node_line)
								output_file.write(sequence_line)
					
				output_file.close()
				fna_file.close()
				gff_file.close()
	#command="rm "+input+"/*.fai"
	#os.system(command)
		
def opts():
	parser = OptionParser()
	#parser.add_option("-u", "--usearch", dest="usearch_path", help="path to USEARCH")
	parser.add_option("-e", "--eggnog", action="store_true", dest="eggnog_run", help ="run eggNOG")
	parser.add_option("-p", "--pilercr", action="store_true",dest="pilercr_run", help ="path to pilercr")
	parser.add_option("-s", "--signalP", action="store_true",dest="signalP_run", help ="path to SignalP")
	#parser.add_option("-hmm", "--hmmtop", action="store_true",dest="hmmtop_run", help ="path to HMMTOP")
	#parser.add_option("-c", "--CARD", dest="card_run", help="path to CARD-rgi")
	parser.add_option("-i", "--input", dest="input_path", help="path to Input")
	return(parser.parse_args())
	
def main():
	temp_dir="../temp_outputs"
	make_temp="mkdir "+temp_dir
	os.system(make_temp)
	options, args = opts()
	eggNOG_run=options.eggnog_run
	pilercr_run=options.pilercr_run
	signalP_run=options.signalP_run
	#hmmtop_run=options.hmmtop_run
	input_path=options.input_path
	#Translate(input_path,temp_dir)
	print(pilercr_run,signalP_run,eggNOG_run)
	
	if pilercr_run != None:
		Pilercr_output_path=Pilercr(input_path,temp_dir)
		formatPilercr(Pilercr_output_path,temp_dir)
	if signalP_run !=None:
		SignalP_output_path=SignalP(input_path,temp_dir)
		formatSignalP("/home/projects/group-c/Team3-FunctionalAnnotation/Outputs/SignalP",temp_dir)
	if eggNOG_run != None:
		output_uclust=uclust(input_path,temp_dir)
		input_eggnog=output_uclust[0]
		input_map=output_uclust[1]
		input_eggnog=format_uclust(input_eggnog,temp_dir)
		cluster_dict=mapNodes(input_map)
		eggNOG_output_path=eggnog(input_eggnog,temp_dir)
		formateggNOG(eggNOG_output_path,cluster_dict,temp_dir)

##Run CARD and VFDB by default 
	CARD_output_path=CARD(input_path,temp_dir)
	
	VFDB_output_path=VFDB(input_path,temp_dir)

##Reformat CARD and VFDB  outputs into gff format to be merged
	formatCARD(CARD_output_path,temp_dir)
	formatVFDB(VFDB_output_path,temp_dir)

##Merge all of the reformated files so there is one full gff file for each sample 
	mergeGff(temp_dir)
	#shutil.rmtree(temp_dir)
	#print(time.time()-start_time)
##Get Fasta files for annotations
	GetFASTA(input_path,"./Outputs/CARD","CARD")
	GetFASTA(input_path,"./Outputs/VFDB","VFDB")
	GetFASTA_all(input_path,"./Outputs/merged")
	shutil.rmtree(temp_dir)
if __name__ == "__main__":
    main()

