
# A - raw data formatting :


*10 species* x *5 individuals* x *multiple runs* (x *multiple fastq files*)

Only the **OBIWAN** runs were used, as **WINDU** runs had troubles during the sequencing. 

Each run has **R1**/**R2**/**R3** files (that can be split in multiple fastq files), where **R1** & **R2** contain paired-end reads and **R3** is identical to **R2** and was therefore not used in the analysis.

All paired-end run use the same insert size : **550 bp**.

Split fastq files corresponding to the same run in an individual were merged together (R1 & R2 being kept apart), leading to a total of :
*10 species* x *5 individuals* x *3 (OBIWAN) runs* x *R1/R2* =  *150 fastq files*

For the correspondance table between original file names from the sequencing facility and standardized file names used in the next steps, see :
[correspondance table (table for download)](./resequencing_samples)


### Species and individual aliases :

* **1_Tdi**   (samples: Tdi_01 / Tdi_02 / Tdi_03 / Tdi_04 / Tdi_05)
* **1_Tps**   (samples: Tps_01 / Tps_02 / Tps_03 / Tps_04 / Tps_05)
* **2_Tcm**   (samples: Tcm_01 / Tcm_02 / Tcm_03 / Tdi_04 / Tcm_05)
* **2_Tsi**   (samples: Tsi_01 / Tce_02 / Tce_03 / Tce_04 / Tce_05)
* **3_Tce**   (samples: Tce_01 / Tdi_02 / Tdi_03 / Tdi_04 / Tdi_05)
* **3_Tms**   (samples: Tms_01 / Tms_02 / Tms_03 / Tms_04 / Tms_05)
* **4_Tbi**   (samples: Tbi_01 / Tbi_02 / Tbi_03 / Tbi_04 / Tbi_05)
* **4_Tte**   (samples: Tte_01 / Tte_02 / Tte_03 / Tte_04 / Tte_05)
* **5_Tge**   (samples: Tge_01 / Tge_02 / Tge_03 / Tge_04 / Tge_05)
* **5_Tpa**   (samples: Tpa_01 / Tpa_02 / Tpa_03 / Tpa_04 / Tpa_05)

