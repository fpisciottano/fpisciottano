 multiPAML2.5 is useful to run multiple PAML analyses given the MSA (.fas) and the corresponding trees (.tre) tagged in the clade to test


 REQUIREMENTS: a MSA file in fasta format for each gene, one or multiple tree files (with .tre extension) for each gene, 'branch-site.ctl' and
                     'branch-site_Ho.ctl' control files and this multiPAML2.5 python script file must be located in the same folder


 EXPECTED OUTPUT: multiPAML2.5 will create the directories structure and manage necesary file to run positive selection tests for all given genes 
                     testing one or many branches/clades for each gene, depending on the tagged trees given for each gene. multiPAML2.5 will crate
                     a directory for each gene and inside it a sub-directory for each branch/clade to be tested for that gene.
                     Output files for each run will be parsed and most relevant results (such as log likelihood values, LRT estimator value, resulting
                     p-value and the table of estimated parametres) for each test will be summarized in a general 'results.txt' output file.
                     Positive results will be additionaly gathered in another ouput file for easy access ('PosSel-tests.txt')


  WARNINGS: this version of multiPAML2 (v2.5) is prepared to run positive selection. If you previously want to run multiple M0 (one-ratio test)
                in order to optimice branch lengthes under GY94, please use multiPAML2.5-M0 script before

                - this script asssumens codeml binary is included in the $PATH and can be invoked from any location

                - MSA files must bear the exact same extension as the one set in the filesExt variable (default: .fas) 

                - in order to be properly assigned to the corresponding gene, tree files must START with the exact same name of the MSA file. After
                     starting with their corresponding MSA file name, tree files should vary in order to show differences between them (e.g., diffent
                     tagged branches or clades)

                - sub-directories of the same gene directory, corresponding to different tested branches/clades for that gene, will be named according
                     to the corresponding tree file name. Specifically they will be named with MSA file name and the last part or the corresponding
                     tree file (everything that appears in the tree file name after the LAST underscore).

                     Recommendation: in order to easyly distinguish among results add the name of the tagged branch/clade at the END of the tree file
                     name AFTER an UNDERSCORE.



multiPAML2.5 was created by Dr. Francisco Pisciottano (fpisciottano@gmail.com) at Biology and Experimental Medicine Institute
- National Scientific and Technical Research Council (IBYME-CONICET), Ciudad Autónoma de Buenos Aires, Argentina