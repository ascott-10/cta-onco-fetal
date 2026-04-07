# python_scripts_common/main.py

from config import *
import sys
from main_process_data import pre_process_project_setup, run_process_pipeline
# main.py
def main():
    if len(sys.argv) < 2:
        print("Usage: main.py process <project_name> | refine <project_name>"); sys.exit(1)

    mode = sys.argv[1]

    try:
        project_name = sys.argv[2]
        paths, project_cfg, global_cfg, samples_to_process = pre_process_project_setup(project_name)
        
        if mode == "process":
            if len(sys.argv) != 3: print("Usage: main.py process <project_name>"); sys.exit(1)
            
            print(f"Running PROCESS mode: {project_name}")
            run_process_pipeline(project_name, paths, project_cfg, global_cfg, samples_to_process)

        elif mode == "refine":
            if len(sys.argv) != 3: print("Usage: main.py refine <project_name>"); sys.exit(1)
           
            print(f"Running REFINE mode: {project_name}")
            #refine_processed_data(project_name, paths, project_cfg, global_cfg, samples_to_process)
        
        else:
            print(f"Unknown mode: {mode}"); sys.exit(1)

    except Exception as e:
        print(f"Pipeline failed with error: {e}")
        import traceback; traceback.print_exc(); sys.exit(1)


if __name__ == "__main__":
    main()



# source myconda; conda activate rnaseq-pipe

# python python_common_scripts/main.py process fetal_gonad
# python python_common_scripts/main.py refine fetal_gonad
# python python_common_scripts/main.py process embryos_mixed
# python python_common_scripts/main.py refine embryos_mixed
# python python_common_scripts/main.py whole_project fetal_gonad
# python python_common_scripts/main.py whole_project embryos_mixed