# Investigating PubMed database
## Introduction
The second project you can choose is an investigation of the scientific literature.
In the data/datasets/NCBI/Pubmed directory on the assemblix computers, you can find a copy of the entire PubMed literature database in XML format.
The project is to use your knowledge of graph theory and parallel processing to investigate the structure of publishing in the scientific world.
Specifically, it answers the following questions:

How large a group of co-authors does the average publication have?
Do authors mostly publish using always the same group of authors?
Do authors mainly reference papers with other authors with whom they've co-authored papers (including themselves)?
What is the distribution in time for citations of papers in general, and for papers with the highest number of citations? Do they differ?
Is there a correlation between citations and the number of keywords that papers share? I.e. papers which share the same subject cite each other more often.
For the most-cited papers (define your own cutoff), is the correlation in shared keywords between them and the papers that cite them different from (5) ?
You should choose or define appropriate measures (statistical, graph theoretical) to answer these questions — and motivate your answer in your report!

Note that it probably helps to parse the actual NCBI XML data only once and then save it in another, smaller format suitable to answering the questions (e.g. a NetworkX datastructure).
