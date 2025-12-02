import os
import sys
import time
import subprocess
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime

class RunfileGenerator:
    def __init__(self, nthreads, dsinit, ffcard):
        with open(dsinit) as f:
            self.dsinit_ = f.readlines()
        with open(ffcard) as f:
            self.ffcard_ = f.readlines()
        self.ffcard_filename_ = ffcard
        self.dsinit_filename_ = dsinit
        self.nthreads_ = nthreads
        self.rung_values_ = dict()

    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_value, traceback):
        if hasattr(self, 'ffcard_names_'):
            print('cleaning up temporary ffcard files...')
            for f in self.ffcard_names_:
                subprocess.run([
                    'rm' , '-rf', f
                ])
        if hasattr(self, 'dsinit_names_'):
            print('cleaning up temporary dsinit files...')
            for f in self.dsinit_names_:
                subprocess.run([
                    'rm' , '-rf', f
                ])
        return False  # propagate exceptions normally
    
    def check_seed_line_(self, line, i):
        parts = line.split()
        if len(parts) == 3 and parts[0] == "RNDM":
            seed1 = int(parts[1])
            seed2 = int(parts[2])
            line = "RNDM %i %i\n"%(seed1+i, seed2+i)
        return line

    def check_ffcard_line_(self, line, ffname):
        parts = line.split()
        if len(parts) == 2 and parts[1][:6] == "FFCARD":
            line = "export FFCARD=\"%s\"\n"%(ffname)
        return line

    def check_rung_line_(self, line, i):
        parts = line.split()
        if len(parts) == 2 and parts[0] == "RUNG":
            rung = int(parts[1])
            rungval = 1000+rung+i-1
            line = "RUNG %i\n"%rungval
            self.rung_values_[i] = rungval
            self.rung_original_ = rung
        return line       

    def write_temp_files(self):
        self.ffcard_names_, self.dsinit_names_ = [], []
        for i in range(1, self.nthreads_+1):
            fn_ffcard = '%.2i_%s' %(i, self.ffcard_filename_)
            fn_dsinit = '%.2i_%s' %(i, self.dsinit_filename_)
            
            self.ffcard_names_.append(fn_ffcard)
            self.dsinit_names_.append(fn_dsinit)
            
            with open(fn_ffcard, 'w') as f:
                for line in self.ffcard_:
                    line1 = self.check_seed_line_(line, i)
                    line2 = self.check_rung_line_(line1,i)
                    f.write(line2)
                    
            with open(fn_dsinit, 'w') as f:
                for line in self.dsinit_:
                    f.write(self.check_ffcard_line_(line, fn_ffcard))

    def get_ffcard_names(self):
        return self.ffcard_names_

    def get_dsinit_names(self):
        return self.dsinit_names_

class RunDsbatch:
    def __init__(self, dsinit):
        self.dsinit_fname_ = dsinit

        self.ffcard_fname_ = None
        with open(self.dsinit_fname_) as f:
            for line in f:
                line = line.strip()
                if line.startswith("export FFCARD="):
                    raw = line.split("=", 1)[1].strip().strip('"')
                    self.ffcard_fname_ = raw.rsplit("/", 1)[-1]
                    break
                
    def run(self, stdout=None, stderr=None):
        cmd = ['./docker_run_dsbatch.sh', self.dsinit_fname_]
        subprocess.run(cmd, stdout=stdout, stderr=stderr)

    def _run_parallel_private(self, nthreads):
        with RunfileGenerator(nthreads, self.dsinit_fname_, self.ffcard_fname_) as rfg:
            rfg.write_temp_files()
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            self.rung_values_ = rfg.rung_values_.copy()
            self.rung_original_ = rfg.rung_original_
        
            def _do_run(i, rfg=rfg, timestamp=timestamp):
                runner_i = RunDsbatch(rfg.dsinit_names_[i])

                self.outdir_ = os.path.join(
                    os.getcwd(), 'dsbatch_parallel_%s'%timestamp
                )
                os.makedirs(self.outdir_, exist_ok=True)
                out_path = '%s/%.2i_dsbatch.out'%(self.outdir_,i)
                with open(out_path, "w") as fout:
                    runner_i.run(stdout=fout,stderr=fout)

            with ThreadPoolExecutor(max_workers=nthreads) as ex:
                return list(ex.map(_do_run, range(nthreads)))

    def run_parallel(self, nthreads):
        print("starting parallel execution (%i threads)...."%nthreads)
        self._run_parallel_private(nthreads)
        for i in self.rung_values_:
            subprocess.run([
                'h2root',
                'dragon%i.hbook'%self.rung_values_[i],
                '%s/dragon%i.root'%(self.outdir_,self.rung_values_[i])
                ],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL
            )
            subprocess.run(
                ['rm', '-f', 'dragon%i.hbook'%self.rung_values_[i]]
            )
        print("done...")
        rung_labels = [self.rung_values_[i] for i in self.rung_values_]        
        print("created ROOT files: %s/dragon%i.root through dragon%i.root"%(
            self.outdir_.split('/')[-1], min(rung_labels), max(rung_labels)))
        print("run this command to merge into a single ROOT file:\n\n"
              "cd %s; hadd merged.root *.root\n"
              "\n(change the output from merged.root to whatever you want it to be),"%(
            self.outdir_.split('/')[-1]))
        # subprocess.run([
        #     'hadd',
        #     '%s/dragon%i.root'%(self.outdir_,self.rung_original_),
        #     '%s/*.root'%(self.outdir_)
        #     ]
        #    #,stdout=subprocess.DEVNULL
        #    #,stderr=subprocess.DEVNULL
        # )
        # for i in self.rung_values_:  
        #     subprocess.run(
        #         ['rm', '-f', '%s/dragon%i.root'%(self.outdir_, self.rung_values_[i])]
        #     )
        # print([
        #     'hadd',
        #     '%s/dragon%i.root'%(self.outdir_,self.rung_original_),
        #     '%s/*.root'%(self.outdir_)
        #     ])


if __name__ == '__main__':
    if(len(sys.argv) < 2):
        print('usage: python %s <dsinit file> [<nthreads=4>]'%sys.argv[0])
    if(len(sys.argv) < 3):
        nthreads = 4
    else:
        nthreads = int(sys.argv[2])
        
    runner = RunDsbatch(sys.argv[1])
    runner.run_parallel(nthreads)
    
          
