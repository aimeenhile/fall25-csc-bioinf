# 1. Repository Setup
- Created a public repository called fall25-csc-bioinf on GitHub with the option not to initialize README or .gitignore
- On personal computer, clone the GitHub repo to Linux subsystem: Set up the repo per instruction: `git clone https://github.com/aimeenhile/fall25-csc-bioinf.git`, then go to the newly created folder: `cd fall25-csc-bioinf` and set up the repo:
```
mkdir -p .github/workflows
mkdir -p week1/{code,data,test}
touch week1/{ai.md,report.md}
```
- Set up the CI: create actions.yml and add the lines provided to it. Then, push the repo to GitHub and make sure all CI passes.
```
cd .github/workflows
touch actions.yml
vim actions.yml
cd fall25-csc-bioinf
git add .
git commit -m "set up CI"
git push origin main
```

# 2. Python Setup & Runs
- Extract the 4 datasets provided into `week1/data`
```
unzip data1.zip -d week1/data/data1
unzip data2.zip -d week1/data/data2
unzip data3.zip -d week1/data/data3
unzip data4.zip -d week1/data/data4
```
- Organize the datasets and the source codes (datasets go to `/data` and `.py` files go to `/code`)
- In directory `week1/code/`, run all the experiements with the source Python code
```
ulimit -s 8192000
python main.py ../data/data1
python main.py ../data/data2
python main.py ../data/data3
python main.py ../data/data4
```

# 3. Condon Conversion & Runs
- Codon version of `main.py` was implemented as `main.codon.py`
```
codon run -release code/main.codon.py ./data
```
