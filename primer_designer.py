import primer3
import re
import csv

default_order = ['gcp', 'gcc', 'temps', 'sizes', 'poly_x'] 
def design(sequence_id, sequence_template, sequence_included_region, insert_site, csv_path, append_to_csv=False, left_search_site=[-700,-100], right_search_site=[100,700], max_tm_diff=2, var_order=default_order): 
    """
    Iteratively attempts to design primers by selectively relaxing set parameters. Results will be printed and exported to a csv designated by csv_path. 
    
    Args: 
        sequence_id: 
            ID information about this DNA sequence
        sequence_template: 
            String containing entire sequence
        sequence_included_region: 
            [start, length] region of sequence to search for primers
        insert_site: 
            index of the insert site to surround with primers
        csv_path: 
            path to csv to add generated primers. will create a new csv unless append_to_csv is set to True. 
        append_to_csv: 
            determines if function creates new csv. if True, please ensure the csv file already exists. default False
        left_search_site: 
            area to search for forward primers in the form [left_index, right_index] such that both indices are relative to the insert site. for example, [-700, 100] would refer to the locus between 700 and 100 bp upstream of the insert site. default [-700, 100]
        right_search_site
            area to search for reverse primers in the form [left_index, right_index] such that both indices are relative to the insert site. for example, [100, 700] would refer to the locus between 100 and 700 bp downstream of the insert site. default [-700, 100]
        max_tm_diff: 
            maximum acceptable difference in melting temperatures between one set of primers. if there is at least one pair of forward/reverse primers that have melting temperatures within this range, the primer set will be accepted. default 2
        var_order:
            the order in which to relax variables, in the form of a list of Strings. default ['gcp', 'gcc', 'poly_x', 'temps', 'sizes'] 
    
    Returns: 
        dictionary containing info for all found primers. if none found, returns None. 
        
    
    """
    
    
    
    # attempt to find primers with default settings
    first_pass = design_primers_specific(sequence_id, sequence_template, sequence_included_region, insert_site, left_search_site=left_search_site, right_search_site=right_search_site)
    left_num_returned = first_pass['PRIMER_LEFT_NUM_RETURNED']
    right_num_returned = first_pass['PRIMER_RIGHT_NUM_RETURNED']
    if left_num_returned > 0 and right_num_returned > 0:
        # check that there is a pair with acceptable melt tm diff
        found_acceptable_pair = False
        for i in range(left_num_returned):
            curr_melt_tm = first_pass[f"PRIMER_LEFT_{i}_TM"]
            for j in range(right_num_returned):
                if abs(first_pass[f"PRIMER_RIGHT_{j}_TM"] - curr_melt_tm) < max_tm_diff:
                    found_acceptable_pair = True
        if found_acceptable_pair: 
            print("Found first pass")
            export_primers("gen_primers.csv", first_pass, append_to_csv)
            return first_pass
        else: 
            print("Found primers, but all had unacceptable melting temp diff.")
        
        
    # default values to be relaxed
    sizes = [22, 23, 25]
    temps = [55, 57, 60]
    gcp = 1.0
    gcc = 2
    poly_x = 2
    
    # define bounds for when program stops. mostly arbitrarily selected, ask ian
    size_lower_lim = 18
    size_upper_lim = 30
    temp_lower_lim = 50
    temp_upper_lim = 65
    gcp_lim = 0.0
    gcc_lim = 0
    poly_x_lim = 8
    
    """
    # used to test end of loop
    size_lower_lim = 22
    size_upper_lim = 25
    temp_lower_lim = 55
    temp_upper_lim = 60
    gcp_lim = 1
    gcc_lim = 2
    poly_x_lim = 2
    """    
    
    # loop through variables until primers are found
    found = False
    nothing_changed = False
    while not found and not nothing_changed:
        nothing_changed = True
        for param in var_order: 
            print(f"ATTEMPTED TO RELAX {param}")
            changed_param = None
            # relax respective variable
            if param == 'gcc': 
                if gcc > gcc_lim:
                    gcc -= 1 
                    changed_param = param
            elif param == 'poly_x': 
                if poly_x < poly_x_lim:
                    poly_x += 1 
                    changed_param = param
            elif param == 'temps':
                if temps[0] > temp_lower_lim:
                    temps[0] -= 1 
                    changed_param = param
                if temps[2] < temp_upper_lim:
                    temps[2] += 1 
                    changed_param = param
            elif param == 'gcp': 
                if gcp > gcp_lim:
                    gcp -= 0.2
                    changed_param = param
            elif param == 'sizes': 
                if sizes[0] > size_lower_lim:
                    sizes[0] -= 1 
                    changed_param = param
                if sizes[2] < size_upper_lim:
                    sizes[2] += 1 
                    changed_param = param
            else: 
                print('loop error')
            print(f"CHANGED {changed_param}")
            if changed_param is not None: 
                nothing_changed = False
            # attempt to find primers with current specific settings
            curr = design_primers_specific(sequence_id, sequence_template, sequence_included_region, insert_site,
                                   left_search_site, right_search_site, sizes=sizes, temps=temps, gcp=gcp, gcc=gcc, poly_x=poly_x)
            print(curr)
            left_num_returned = curr['PRIMER_LEFT_NUM_RETURNED']
            right_num_returned = curr['PRIMER_RIGHT_NUM_RETURNED']
            if left_num_returned > 0 and right_num_returned > 0:
                # check that there is a pair with acceptable melt tm diff
                found_acceptable_pair = False
                for i in range(left_num_returned):
                    curr_melt_tm = curr[f"PRIMER_LEFT_{i}_TM"]
                    for j in range(right_num_returned):
                        if abs(curr[f"PRIMER_RIGHT_{j}_TM"] - curr_melt_tm) < max_tm_diff:
                            found_acceptable_pair = True
                if found_acceptable_pair: 
                    found = True
                    print("FOUND")
                    export_primers(csv_path, curr, append_to_csv)
                    return curr
                else: 
                    print("Found primers, but all had unacceptable melting temp diff.")
    print("COULD NOT FIND SUITABLE PRIMERS")
    
    
    
    
def design_primers_specific(sequence_id, sequence_template, sequence_included_region, 
                            insert_site, left_search_site=[-700,-100], right_search_site=[100,700],
                            sizes=[22, 23, 25], temps=[55,57,60], gcp=1.0, gcc=2, poly_x=2):
    """
    Helper function for iterative_design(). Attempts to design primers with specific parameter settings using Primer3. 
    """
    sizes.sort()
    temps.sort()
    primer_dict = primer3.bindings.designPrimers(
        {
            'SEQUENCE_ID': sequence_id,
            'SEQUENCE_TEMPLATE': sequence_template,
            'SEQUENCE_INCLUDED_REGION': sequence_included_region,
            'SEQUENCE_PRIMER_PAIR_OK_REGION_LIST': [[insert_site+left_search_site[0], 
                                                     abs(left_search_site[0]-left_search_site[1]), 
                                                     insert_site+right_search_site[0], 
                                                     abs(right_search_site[1]-right_search_site[0])]]
        }, 
        {
            'PRIMER_TASK': 'pick_primer_list',
            'PRIMER_NUM_RETURN': 30,
            'PRIMER_MIN_SIZE': sizes[0],
            'PRIMER_OPT_SIZE': sizes[1],
            'PRIMER_MAX_SIZE': sizes[2],
            'PRIMER_MIN_TM': temps[0],
            'PRIMER_OPT_TM': temps[1],
            'PRIMER_MAX_TM': temps[2],
            'PRIMER_OPT_GC_PERCENT': 50.0,
            'PRIMER_WT_GC_PERCENT_LT': gcp,
            'PRIMER_WT_GC_PERCENT_GT': gcp,
            'PRIMER_GC_CLAMP': gcc,
            'PRIMER_MAX_POLY_X': poly_x
        })
    return primer_dict



def export_primers(output_path, primer_dict, append): 
    """
    Helper function for iterative_design(). Exports primers to specified csv in output_path.
    """
    
    csv_file = open(output_path, 'a' if append else 'w')
    csvwriter = csv.writer(csv_file)
    csvwriter.writerow(['name', 'sequence'])
    
    primers_matrix = []
    for key in primer_dict:
        if re.fullmatch(r"PRIMER_((LEFT)|(RIGHT))_[0-9]+_SEQUENCE", key):
            primers_matrix.append([key[:-9], primer_dict[key]])
    
    csvwriter.writerows(primers_matrix)
    csv_file.close()
    
    
    
def get_oligos_specific_regions(sequence_id, sequence_template, csv_path, left_search_site, right_search_site, append_to_csv=False, max_tm_diff=2, var_order=default_order): 
    
    
    left_section = sequence_template[left_search_site[0]:left_search_site[1]]
    right_section = sequence_template[right_search_site[0]:right_search_site[1]]
    
    middle_cut_out = left_section+right_section
    
    print(design(sequence_id, middle_cut_out, (0, len(middle_cut_out)), 0, csv_path, append_to_csv=append_to_csv, left_search_site=left_search_site,
           
           right_search_site=[len(middle_cut_out)-(right_search_site[1]-right_search_site[0]),len(middle_cut_out)], 
           
           max_tm_diff=max_tm_diff, var_order=var_order))
    
    """left_results = design(sequence_id, left_section, (0, len(left_section)), 0, csv_path, append_to_csv=append_to_csv, left_search_site = left_search_site, right_search_site=[0,0], max_tm_diff=max_tm_diff, var_order=var_order)
    
    right_results = design(sequence_id, right_section, (0, len(right_section)), 0, csv_path, append_to_csv=append_to_csv, left_search_site = [0,0], right_search_site=[0,right_search_site[1]-right_search_site[0]], max_tm_diff=max_tm_diff, var_order=var_order)"""
    
    

    
    
def prune_results(orig_csv_path, new_csv_path, left_indices, right_indices, append=True): 
    """
    Prunes results of primers in orig_csv_path by putting them in another csv designated by new_csv_path. Meant to be used after visual validation by the user to determine which generated primers to keep. 
    
    Args: 
        orig_csv_path: 
            csv with all primers to choose from. typically will be the generated csv from iterative_design()
        new_csv_path: 
            csv to put selected primers in. if append is False, the function will create a new csv file, overwriting the original contents. 
        left_indices: 
            list of forward primers to keep
        right_indices: 
            list of reverse primers to keep
        append
            set to True if you want to append primers to an existing csv. default False (will create new csv file or overwrite old contents)
    """
    
    orig_csv_file = open(orig_csv_path)
    csvreader = csv.reader(orig_csv_file)
    next(csvreader)
    
    new_csv_file = open(new_csv_path, 'a' if append else 'w')
    csvwriter = csv.writer(new_csv_file)
    if not append: 
        csvwriter.writerow(['name', 'sequence'])
    all_primers = {}
    for row in csvreader: 
        all_primers[row[0]] = row[1]
    orig_csv_file.close()
    new_matrix = []
    for left_index in left_indices: 
        
        name = f"PRIMER_LEFT_{left_index}"
        
        if name not in all_primers.keys(): 
            raise Exception('primer not in list')
            
        new_matrix.append([name, all_primers[name]])
        
    for right_index in right_indices: 
        name = f"PRIMER_RIGHT_{right_index}"
        
        if name not in all_primers.keys(): 
            raise Exception('primer not in list')
        
        new_matrix.append([name, all_primers[name]])
    
    csvwriter.writerows(new_matrix)
    new_csv_file.close()
    
    