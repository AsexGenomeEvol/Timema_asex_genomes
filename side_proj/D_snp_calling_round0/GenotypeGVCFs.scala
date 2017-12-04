package test.ug

import org.broadinstitute.gatk.queue.extensions.gatk.GenotypeGVCFs
import org.broadinstitute.gatk.queue.QScript
import org.broadinstitute.gatk.queue.extensions.gatk._
/**
* An example building on the intro ExampleCountReads.scala.
* Runs an INCOMPLETE variant calling pipeline with just the UnifiedGenotyper, VariantEval and optional VariantFiltration.
* For a complete description of the suggested for a variant calling pipeline see the latest version of the Best Practice Variant Detection document
*/

class ExampleUnifiedGenotyper extends QScript {
// Create an alias 'qscript' to be able to access variables
// in the ExampleUnifiedGenotyper.
// 'qscript' is now the same as 'ExampleUnifiedGenotyper.this'

qscript =>
// Required arguments. All initialized to empty values.

@Input(doc="The reference file for the bam files.", shortName="R")
var referenceFile: File = _ // _ is scala shorthand for null

@Input(doc="Gvcf file list to genotype.", shortName="V")
var gvcfsFile: File = _
// The following arguments are all optional.

@Input(doc="Output filename.", shortName="o")
var outFile: File = _

@Input(doc="An optional file with a list of intervals to proccess.", shortName="L", required=false)
var intervals: File = _

@Input(doc="An optional DBSNP File.", shortName="D", required=false)
var dbsnps: File = _

// This trait allows us set the variables below in one place,
// and then reuse this trait on each CommandLineGATK function below.

trait UnifiedGenotyperArguments extends CommandLineGATK {
this.reference_sequence = qscript.referenceFile
this.intervals = if (qscript.intervals == null) Nil else List(qscript.intervals)

// Set the memory limit to 2 gigabytes on each command.
this.memoryLimit = 40
}

def script() {

// Create the four functions that we may run depending on options.
val genotyper = new GenotypeGVCFs with UnifiedGenotyperArguments

genotyper.scatterCount = 40
genotyper.variant = Seq(qscript.gvcfsFile)
genotyper.dbsnp = qscript.dbsnps
genotyper.out = qscript.outFile
genotyper.nt = 1
add(genotyper)
}
}

