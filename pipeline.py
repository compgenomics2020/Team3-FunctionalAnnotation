#!/usr/bin/env python2
from optparse import OptionParser
import os
from os import popen
import shutil
import subprocess
import multiprocessing

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
def CARD(input_dir,output_dir):
	print("Running CARD-rgi......")
	make_temp="mkdir "+output_dir+"/CARD"
	os.system(make_temp)
	CARD_output_path=output_dir+"/CARD"
	#command="ls "+input_dir+" -p . | grep -v /$"

	#filename = popen(command,"r").readlines()
	#filename = filename[19:]
	#filename1 = []
	#for name in filename:
		#filename1.append(name[:-1])
	#print(filename1)
	counter = 0
	#for files in filename1:
	for files in os.listdir(input_dir):
		prefix=files.split("_")[0]
		#prefix = files[0:7]
		command = "rgi -i "+input_dir+"/"+files+" -o "+CARD_output_path+"/"+prefix+"_CARD-RGI_coding"
		subprocess.call(command,shell=True)
		counter +=1
		print("Done ...",counter)
	return(CARD_output_path)

def VFDB(input_dir,output_dir):
	print("Running blastp on VFDB......")
	make_temp="mkdir "+output_dir+"/VFDB"
	os.system(make_temp)
	VFDB_output_path=output_dir+"/VFDB"
	#filename = popen("ls ../gene_prediction_updated/Coding_FASTA -p . | grep -v /$","r").readlines()
	#filename = filename[15:]
	#filename1 = []
	#for name in filename:
		#filename1.append(name[:-1])

	#for files in filename1:
	for files in os.listdir(input_dir):
		prefix=files.split("_")[0]
		#prefix = files[0:7]
		command = "blastn -db ../Tools/VFDB/Virulence_Factors_core -query "+input_dir+"/"+files+" -out "+VFDB_output_path+"/"+prefix+"_VFDB_coding -max_hsps 1 -max_target_seqs 1 -num_threads 4 -evalue 1e-5 "
		subprocess.call(command,shell=True)
	return(VFDB_output_path)
def Pilercr(input_dir,output_dir)
	make_temp="mkdir "+output_dir+"/pilercr"
	os.system(make_temp)
	pilercr_output_path=output_dir+"/pilercr"
	for filename in os.listdir(input_dir):
		prefix = filename.split("_")[0]
		file_path=input_dir+"/pilercr/"+filename
		command = "../Tools/pilercr/pilercr1.06/pilercr -in ../gene_prediction_updated/Coding_FASTA/"+file+" -out "+pilercr_output_path+"/"+prefix+"_pilercr_coding"
		subprocess.call(command,shell=True)
def eggnog(input,output_dir):
	make_temp="mkdir "+output_dir+"/eggNOG"
        os.system(make_temp)
	eggNOG_annotation_path=output_dir+"/eggNOG"
	call_eggNog="emapper.py -i "+input+"  --data_dir ../Tools/eggnog_db -m diamond --dmnd_db /home/projects/group-c/Team3-FunctionalAnnotation/Tools/eggnog_db/eggnog_proteins.dmnd --translate -d bact --output_dir "+eggNOG_annotation_path+" -o pipeline"
	print(call_eggNog)
	os.system(call_eggNog)
	return(eggNOG_annotation_path)
def SignalP(input_dir,output_dir):
	make_temp="mkdir "+output_dir+"/SignalP"
	os.system(make_temp)
	SignalP_path=output_dir+"/signalP"
	for filename in os.listdir(input_dir):
		sample=filename.split("_")[0]
    		os.system("/home/projects/group-c/Team3-FunctionalAnnotation/Tools/signalp-5.0b/bin/signalp -fasta "+input_dir+"/"+filename+" -org gram+ -format short -gff3 -prefix "+output_dir+"/SignalP/"+sample+"_signalp")
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
def formatCARD(input_dir,output_dir):
	make_temp="mkdir -p "+output_dir+"/format/CARD"
	os.system(make_temp)
	path_to_CARD=input_dir
	print("formatting...card")
	for filename in os.listdir(path_to_CARD):
		#print("here1")
		if "gff" in filename:
			#print("here2")
			sample=filename.split("_")[0]
			path_CARD_file=path_to_CARD+"/"+filename
			CARD_file=open(path_CARD_file,"r")
			CARD_output_path_gff=output_dir+"/format/CARD/"+sample+"_card_mapped.gff"
			CARD_output_gff=open(CARD_output_path_gff,"w+")
			for line in CARD_file:
				if line.startswith("NODE"):
					#print("here3")
					annotation=line.split("\t")
					node=annotation[0]
					node=node.split("_")
					node=node[:-1]
					node_join="_"
					node=node_join.join(node)
					#print(node)
					annotation[0]=node
					annotation_join="\t"
					print(annotation)
					annotation=annotation_join.join(annotation)
					gff_line=annotation
					#print(gff_line)
					CARD_output_gff.write(gff_line)
			CARD_output_gff.close()	
	
def formatVFDB(input_dir,output_dir):
	VFDB_dir=input_dir
	make_temp="mkdir -p "+output_dir+"/format/VFDB"
	os.system(make_temp)

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
def formatSignalP(input_dir,output_dir):
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
def mergeGff(output_dir):
	make_temp="mkdir -p ./Outputs/merged"
	os.system(make_temp)
	sample_dict={}
	mapped_dir_path=output_dir+"/format"
	for dir_name in os.listdir(mapped_dir_path):
		if dir_name=="eggNog" or dir_name=="CARD" or dir_name=="VFDB" or dir_name=="eggNOG" or dir_name=="SignalP":
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
			print(item)
			value_sort=sorted(item, key=lambda key1: (int((key1.split(":")[1]).split("-")[1])))
			for get_item in value_sort:
				loop_lines=feature_dict.get(get_item)
				node=get_item
				header="##sequence_region "+node+"\n"
				output_file.write(header)
				for annotation in loop_lines:
					output_file.write(annotation)
		output_file.close()
	
def opts():
	parser = OptionParser()
	parser.add_option("-u", "--usearch", dest="usearch_path", help="path to USEARCH")
	#parser.add_option("-e", "--eggnog", dest="eggnog_path", help ="path to eggNog")
	#parser.add_option("-c", "--CARD", dest="card_path", help="path to CARD-rgi")
	parser.add_option("-i", "--input", dest="input_path", help="path to Input")
	return(parser.parse_args())
	
def main():
	temp_dir="../temp_outputs"
	make_temp="mkdir "+temp_dir
	os.system(make_temp)
	options, args = opts()
	input_path=options.input_path
	Translate(input_path,temp_dir)
	#output_uclust=uclust(input_path,temp_dir)
	#input_eggnog=output_uclust[0]
	#input_map=output_uclust[1]
	process_CARD=multiprocessing.Process(target=CARD, args=(input_path,temp_dir,))
	CARD_output_path=output_dir+"/CARD"
	process_VFDB=multiprocessing.Process(target=VFDB, args=(input_path,temp_dir,))
	VFDB_output_path=output_dir+"/VFDB"
	process_CARD.start()
	process_VFDB.start()
	process_CARD.join()
	process_VFDB.join()
	#SignalP_output_path=SignalP(input_path,temp_dir)
	#eggnog(input_eggnog,temp_dir)
	#cluster_dict=mapNodes(input_map)
	#formateggNOG("../Outputs/eggNOG/Test_annotations.emapper.annotations",cluster_dict,temp_dir)
	#formatCARD(CARD_output_path,temp_dir)
	#formatVFDB(VFDB_output_path,temp_dir)
	#formatSignalP("/home/projects/group-c/Team3-FunctionalAnnotation/Outputs/SignalP",temp_dir)
	#mergeGff(temp_dir)
	#shutil.rmtree(temp_dir)
if __name__ == "__main__":
    main()

