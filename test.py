if len(good_drafts.keys()) > 0:
        
    if __name__ == '__main__':  
        Parallel(n_jobs = -1, verbose = 5)(delayed(search_16S)
        (ref_dir + 'user/' + domain, d) for d in good_drafts)

with open(ref_dir + 'user/' + domain + '/' + 'draft.combined_16S.fasta', 'w') as draft_fasta_out:        
        for d in good_drafts.keys():
            
            count_16S = subprocess.Popen('grep -c \'>\' ' + ref_dir + 'user/' + domain + '/' + d + '/' + d + '.16S.fasta', shell = True, executable = executable, stdout = subprocess.PIPE)
            n16S = count_16S.communicate()[0]
            n16S = n16S.rstrip()
            n16S = int(n16S)
            
            ## Were any 16S rRNA genes found in the draft assembly?
            
            if n16S > 0:
                
                ## Make an entry for the draft in summary_complete
                
                summary_complete.loc[d, 'organism_name'] = d + '_' + 'DRAFT' + '_' + good_drafts[d]
                summary_complete.loc[d, 'assembly_level'] = 'Draft'
                
                ## Copy the directory to refseq so that pathway-tools can find it.
                
                cp = subprocess.Popen('cp -r ' + ref_dir + 'user/' + domain + '/' + d + ' ' + ref_dir_domain + 'refseq/', shell = True, executable = executable)
                
                ## Print one of the 16S rRNA genes to combined_16.
                
                keep = True
                                    
                for record in SeqIO.parse(ref_dir + 'user/' + domain + '/' + d + '/' + d + '.16S.fasta', 'fasta'):
                    if keep == True:
                        
                        new_record = SeqRecord(record.seq)                        
                        new_record.id = good_drafts[d]                       
                        new_record.description = ''
                        
                        SeqIO.write(new_record, draft_fasta_out, 'fasta')                        
                        keep = False
                        
                    else:
                        continue
                    
            else:
                print 'sorry, no 16S rRNA genes found in draft assembly', d
