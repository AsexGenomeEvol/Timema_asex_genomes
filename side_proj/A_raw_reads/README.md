
# A - raw data formatting :

(go to next step : [read trimming](../B_cleaned_reads))

-----------

*10 species* x *5 individuals* x *multiple runs* (x *multiple fastq files*)

Only the **OBIWAN** runs were used, as **WINDU** runs had troubles during the sequencing (only its **R1** reads are reliable). 

Each run has **R1**/**R2**/**R3** files (that can be split in multiple fastq files), where **R1** & **R2** contain paired-end reads and **R3** is identical to **R2** and was therefore not used in the analysis.

All paired-end runs use the same insert size : **550 bp** with reads of length : **100 bp**.

Split fastq files corresponding to the same run of an individual were merged together (R1 & R2 being kept apart), leading to a total of :
*10 species* x *5 individuals* x *3 (OBIWAN) runs* x *R1/R2* =  *150 fastq files*

For the correspondance table between original file names from the sequencing facility and standardized file names used in the next steps, see :
[correspondance table (table for download)](./resequencing_samples)

The paired-end fastq files were then trimmed for presence of Illumina adapters and quality [here](../B_cleaned_reads).


### Species and individual aliases :

* **1_Tdi**   (samples: *Tdi_01* to *Tdi_05*)
* **1_Tps**   (samples: *Tps_01* to *Tps_05*)
* **2_Tcm**   (samples: *Tcm_01* to *Tcm_05*)
* **2_Tsi**   (samples: *Tsi_01* to *Tce_05*)
* **3_Tce**   (samples: *Tce_01* to *Tdi_05*)
* **3_Tms**   (samples: *Tms_01* to *Tms_05*)
* **4_Tbi**   (samples: *Tbi_01* to *Tbi_05*)
* **4_Tte**   (samples: *Tte_01* to *Tte_05*)
* **5_Tge**   (samples: *Tge_01* to *Tge_05*)
* **5_Tpa**   (samples: *Tpa_01* to *Tpa_05*)


```
               ___ T. poppensis      1_Tps ♀♂  = sexual
              /  
             /\___ T. douglasi       1_Tdi ♀   = asexual
            /
           /\  ___ T. californicum   2_Tcm ♀♂
          /  \/
         /    \___ T. shepardi       2_Tsi ♀
        /\
       /  \   ____ T. cristinae      3_Tce ♀♂
      /    \ /
      \     \_____ T. monikensis     3_Tms ♀
       \
        \      ___ T. bartmani       4_Tbi ♀♂
         \    /
          \  /\___ T. tahoe          4_Tte ♀
           \/
            \  ___ T. podura         5_Tpa ♀♂
             \/  
              \___ T. genevieve      5_Tge ♀
              
(from Kamil page: https://github.com/AsexGenomeEvol/timema_assembly)
```

