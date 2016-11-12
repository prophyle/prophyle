# How to run experiments

This procedure should be performed whenever a major
update is going to be merged to master.

1. Clean all subdirectories
   
  ```
  make clean
  ```

2. Adjust parts of Makefile responsible for logs if it is necessary (e.g., in case of a new feature)

3. Activate the Conda environment

  ```
  conda activate prophyle
  ```
  
4. Run the experiments (with the same number as your machine has)

  ```
  make -j 24
  ```

5. Collect logs
  ```
  ./view_logs_md.sh > _logs/paprika_NB_BRANCHENAME.md
  ./view_logs_txt.sh > _logs/paprika_NB_BRANCHENAME.txt
  ```

6. Create a git commit and push it to the repository

7. Verify that RPM and memory footprint are OK

8. Merge with master
