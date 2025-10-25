# Week 4 Deliverable

Time needed to complete: 24-30 hours

1. Set up
- Collect the data from [ksw2 directory](https://github.com/lh3/ksw2/blob/master)
- Put data files into the `week4/data` folder
- Set up the Bash script `evaluate.sh`
- Set up `utils.py` to read input from FASTA files 
- Set up the CI `actions.yml` file using the format from [https://github.com/inumanag/fall25-csc-bioinf/blob/main/.github/workflows/actions.yml](https://github.com/inumanag/fall25-csc-bioinf/blob/main/.github/workflows/actions.yml)
- Setting up `python/utils.py` to read data from the FASTA files
- Setting up the `python/main.py` so that it will print the output to the CI to test. 

2. Review methods
- Review lecture notes to implement the four algorithms
- I also looked at the resources from [Langmead Lab](https://www.langmead-lab.org/teaching.html) to implement the methods:
    - [global alignment](https://www.cs.jhu.edu/~langmea/resources/lecture_notes/12_global_alignment_v2.pdf)
    - [global alignment](https://github.com/BenLangmead/ads1-notebooks/blob/master/3.02_GlobalAlignment.ipynb)
    - [local alignment](https://www.cs.jhu.edu/~langmea/resources/lecture_notes/13_local_alignment_v2.pdf)
    - [edit distance](https://github.com/BenLangmead/ads1-notebooks/blob/master/3.01_EditDistanceDP.ipynb)
    - [semi-global](https://www.cs.cmu.edu/~durand/03-711/2023/Lectures/20231001_semi-global.pdf)

3. Implement Python codes
- I first implemented the global, local, and affine alignment, then lastly, the semi-global alignment 
- At first I stored the distance and trackback in 2 separate full matrix so the runtime for Python methods was extremely large. I had to ask ChatGPT to suggest a way to efficiently implement the complete traceback.
- I also ran into errors when implementing the semi-global alignment method as I was confused on how to implement it and I was missing a lot of cases for that method, so I looked up online and ask AI
- The score for the semi-global method for q5 vs. t5 was negative so I also asked ChatGPT why it was so low and tried to fix it, but the source file also have the t5 marked as "trigger bug" so I ignore it after many tries. 
- Initally, the Python code took almost 30 minutes to run each time so I keep going back and forth on how to optimize the runtime.
- I gave up trying when the runtime for all four Python methods was 20 minutes.

4. Implement Codon code
- Modified `python/utils.py` and `python/main.py`to Codon. I ran into PyError during the process so I also had to go back and forth to modify the files and fix the errors, but then I found out it was because my careless mistake that caused the error, but eventually it got fixed.
- Modified the Python code to Codon code (many trials and errors here)
- Testing the Codon code and compare the results for q1-q5 vs. t1-t5 for Python and Codon 