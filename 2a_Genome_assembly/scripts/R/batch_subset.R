batch_subset <- function(scf_asm, batch){
  if(batch == 1){
    scf_asm_batch <- scf_asm[scf_asm$fasteris == F &
                             scf_asm$fewdata == T &
                             scf_asm$mse == F &
                             scf_asm$pse == T &
                             scf_asm$mpe == F &
                             scf_asm$soft == 'SOAP', ]
  } else if(batch == 2){
    scf_asm_batch <- scf_asm[scf_asm$fasteris == F &
                             scf_asm$fewdata == F &
                             scf_asm$mse == F &
                             scf_asm$pse == T &
                             scf_asm$mpe == F &
                             scf_asm$soft == 'SOAP', ]
  } else if(batch == 3){
    scf_asm_batch <- scf_asm[scf_asm$fasteris == F &
                             scf_asm$fewdata == T &
                             scf_asm$mse == F &
                             scf_asm$pse == T &
                             scf_asm$mpe == F &
                             c(scf_asm$soft == 'abyss' | scf_asm$soft == 'BESST'), ]
  } else {
    scf_asm_batch <- scf_asm[NULL, ]
  }
  return(scf_asm_batch)
}
