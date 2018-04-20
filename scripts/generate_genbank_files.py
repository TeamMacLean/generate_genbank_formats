#generate genbank files, one for each row, in an excel file.

import sys, os, time

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.SeqFeature import SeqFeature, FeatureLocation



with open(sys.argv[1]) as inputfilehandler:
    #read the header out
    inputfilehandler.readline()

    for line in inputfilehandler:
        #read each line/row at a time, each row has all required info to create a .gbk file
        line=line.rstrip()
        if line == "":
            continue
        else:
            array=line.split("\t")
            EC_assembly = array[1]
            backbone = array[2]
            backbone_sequence = array[4]
            backbone_start = array[5]
            backbone_end = array[6]

            EC1 = array[7]
            EC2 = array[8]
            EC3 = array[9]

            EC1_seq = array[10]
            EC1_digested = array[14]
            EC2_seq = array[15]
            EC2_digested = array[19]
            EC2_digested_start = array[20]
            EC2_digested_end = array[21]

            EC3_seq = array[22]
            EC3_digested = array[22]
            EC3_digested_start = array[27]
            EC3_digested_end = array[28]
            final_sequence = array[29].lower()
            EC1_start_position = array[32]
            EC1_end_position = array[33]
            EC2_start_position = array[34]
            EC2_end_position = array[35]
            EC3_start_position = array[36]
            EC3_end_position = array[37]
            final_sequence_overlap_start = array[39]
            final_sequence_overlap_end = array[40]



            #print EC_assembly, backbone, backbone_sequence, EC1, EC2, EC3, EC1_seq, EC1_digested,EC2_seq, EC2_digested, EC3_seq, EC3_digested, final_sequence, EC1_start_position, EC1_end_position, EC2_start_position, EC2_end_position, EC3_start_position, EC3_end_position
            #print EC_assembly, EC1_start_position, EC1_end_position, EC2_start_position, EC2_end_position, EC3_start_position, EC3_end_position

            # Create a sequence
            sequence_object = Seq(final_sequence, IUPAC.unambiguous_dna)

            # Create a record
            record = SeqRecord(sequence_object,
                               id=EC_assembly + ".1", # random accession number
                               name=EC_assembly)

            backbone_feature = SeqFeature(FeatureLocation(0, len(backbone_sequence)), type="background")
            backbone_feature.qualifiers["label"] = "vector - " + backbone
            record.features.append(backbone_feature)

            seq_description = None
            if not EC1 == "" and EC1_start_position != "":
                seq_description = "EC used in the circular construct; EC1 - " + EC1
                EC1feature=SeqFeature(FeatureLocation(start=int(EC1_start_position.strip())-1, end=int(EC1_end_position.strip())), type="CDS")
                EC1feature.qualifiers["label"] = EC1
                record.features.append(EC1feature)

            if not EC2 == "" and EC2_start_position != "":
                seq_description += ", EC2 - " + EC2
                EC2feature=SeqFeature(FeatureLocation(start=int(EC2_start_position.strip())-1, end=int(EC2_end_position)), type="CDS")
                EC2feature.qualifiers["label"] = EC2
                record.features.append(EC2feature)
            if not EC3 == "" and EC3_start_position != "":
                seq_description += ", EC3 - " + EC3
                EC3feature=SeqFeature(FeatureLocation(start=int(EC3_start_position)-1, end=int(EC3_end_position)), type="CDS")
                EC3feature.qualifiers["label"] = EC3
                record.features.append(EC3feature)


            record.description=seq_description
            date_today = str(time.strftime("%d-%b-%Y").strip()).upper()
            record.annotations ={"organism": "Cloning Vector", "data_file_division":"SYN", "date":date_today, "residue_type": "DNA"}
            base_count=str(final_sequence.count("a")) + " a      \t"+ str(final_sequence.count("c")) + " c      \t" + str(final_sequence.count("g")) + " g      \t" + str(final_sequence.count("t")) + " t"

            record.annotations["source"] = "vector"
            record.annotations["gi"] = EC_assembly
            record.annotations["base count"] = base_count

            if array[5] != "" and array[6] != "":
                misc = SeqFeature(FeatureLocation(int(array[5])-1, int(array[6])), type="misc_feature")
                misc.qualifiers["label"]="4bp\overhang"
                record.features.append(misc)
            if array[20] != "" and array[21] != "":
                misc = SeqFeature(FeatureLocation(int(array[20])-1, int(array[21])), type="misc_feature")
                misc.qualifiers["label"] = "4bp\overlap"
                record.features.append(misc)
            if array[68] != "" and array[69] != "":
                misc = SeqFeature(FeatureLocation(int(array[68])-1, int(array[69])), type="misc_feature")
                misc.qualifiers["label"] = "3xFlag"
                record.features.append(misc)
            if array[70] != "" and array[71] != "":
                misc = SeqFeature(FeatureLocation(int(array[70])-1, int(array[71])), type="misc_feature")
                misc.qualifiers["label"] = "rbcS-Ter"
                record.features.append(misc)
            if array[72] != "" and array[73] != "":
                misc = SeqFeature(FeatureLocation(int(array[72])-1, int(array[73])), type="misc_feature")
                misc.qualifiers["label"] = "amp/carben_R"
                record.features.append(misc)
            if array[74] != "" and array[75] != "":
                misc = SeqFeature(FeatureLocation(int(array[74])-1, int(array[75])), type="misc_feature")
                misc.qualifiers["label"] = "2x35S_promoter"
                record.features.append(misc)
            if array[76] != "" and array[77] != "":
                misc = SeqFeature(FeatureLocation(int(array[76])-1, int(array[77])), type="misc_feature")
                misc.qualifiers["label"] = "TMV_Omega"
                record.features.append(misc)
                '''
            if array[27] != "" and array[28] != "":
                misc = SeqFeature(FeatureLocation(int(array[27])-1, int(array[28])), type="misc_feature")
                misc.qualifiers["label"] = "4bp\overlap"
                record.features.append(misc)
            if array[39] != "" and array[40] != "":
                misc = SeqFeature(FeatureLocation(int(array[39])-1, int(array[40])), type="misc_feature")
                misc.qualifiers["label"] = "4bp\overlap"
                record.features.append(misc)
                '''

            #record.annotations.keys()  are ['comment', 'sequence_version', 'source', 'taxonomy', 'keywords', 'references','accessions', 'data_file_division', 'date', 'organism', 'gi']
            # Save as GenBank file

            output_file = open(EC_assembly + '.gb', 'w')
            SeqIO.write(record, output_file, 'genbank')
            output_file.close()
            cmd="bash scripts/add_circular_basecounts.sh " + EC_assembly + ".gb \"" + base_count + "\" > " + "temp; mv temp " + EC_assembly + ".gb"
            print(cmd)
            os.system(cmd)


exit(0)
