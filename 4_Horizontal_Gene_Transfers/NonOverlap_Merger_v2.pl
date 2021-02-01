# courtesy of Paul simion !!!!!!!!!!!!!!!!!!!!!!!
# 27/10/2014
# Objectif: Ouvrir un alignement, comparer les séquences partielles d'une même espèce, vérifier leur overlap et concaténer/supprimer selon les cas.			
# REQUIREMENTS: need sequences sur une seule ligne / meme ID (au moins avant le premier '| ')
#  ATTENTION: il faut un pipe dans le nom des séqs !
#  ATTENTION: for now, les seuils pour les chevauchements sont: si longueur chevauchement > 25 (can be modified below; peut aussi incorporer le nbre de diffs max)

# Usage: 3 arguments:
# perl NonOverlap_Merger_v2.pl nom_alignement outfile_basename notification
# exemple: perl ~/scripts/NonOverlap_Merger_v2.pl infile.fasta outfile notification
# le(s) seq(s) consensus sont output dans 'outfile.out'

#! /usr/bin/perl

use strict;
use warnings;
#use lib '/home/paul/0_Paul/Programmation/scripts_perl_Paul';
use List::MoreUtils ':all';

### déclaration des variables
my @ali = ();
my @matrix = ();
my @species_names = ();
my @seq_names = ();
my @sequences = ();
my @keys = ();
my @sptoconcatenate = ();
my @tmp_consensus = ();
my @tmp_present = ();
my $line = '';
my $current_contig = '';
my $nb_seq = '';
my $max_length = '';
my $seq = '';
my $flag = '';
my $best_contig = '';
my %bilan = ();
my %seqtokeep = ();
my %seqtoconcatenate = ();
my %sptoconcatenate = ();
my %seq_consensus = ();
my %longueur_seq = ();
my %sequences = ();
my %difference = ();
my %chevauchement = ();
my %identity = ();
my @ctgs = ();
my %alerte = ();
my %list_contigs = ();
my $percentid = '';

### vérification ouverture script et information utilisateur
print "Alignment name : $ARGV[0]\nOutfiles basename : $ARGV[1]\n"; 

### ouverture alignement
open (INFILE0,$ARGV[0]) or die "\nImpossible d'ouvrir le fichier $ARGV[0]";
@ali = <INFILE0>;
close (INFILE0);

### lecture des noms des séquences ET des espèces => dans deux tableaux
# (ATTENTION au délimiteur de champs !!!!)
foreach my $line (@ali) {
	if ( $line =~ /^>([^\|]*)\|.*$/ ) {
		push (@species_names, $1);			# on garde le nom de l'espèce
		$line =~ s/\n//;
		push (@seq_names, $line);			# on garde le nom du contig
		#print"\nofficiel: $line";
	}
}

### lecture des séquences => dans une matrice
foreach my $line (@ali) {
	if ( $line =~ /^>([^\|]*)\|.*$/ ) {
		$line =~ s/\n//;
		#$line =~ s/>//;
		$current_contig = $line;
		#print "\n\n##$current_contig";
	}
	elsif ( $line =~ /^[^>]/ ) {
 		my @splited_line = split ('', $line); 			# on morcelle la ligne
 		push (@matrix, [@splited_line]);				# on continue la matrice générale
		$sequences{$current_contig}=$line;				# on sauvegarde le couple contig => sequence
		#$longueur_raw_seq{$current_contig}=$#splited_line;	# on sauvegarde le couple contig => longueur sequence alignée
		@tmp_present=();
		for (my $x = 0; $x <= "$#splited_line"; $x++) {
			if ($splited_line[$x] !~ m/[-\?X\*]/g) {push( @tmp_present ,$splited_line[$x]);	}}	# on calcule le nombre de caractères renseignés
		$longueur_seq{$current_contig}=scalar(@tmp_present);	# on sauvegarde le couple contig => longueur sequence renseignée
		print "\nlecture de $current_contig ($longueur_seq{$current_contig} AA)";
	}
}

### mesures des dimensions de la matrice
$nb_seq = $#matrix+1;
for (my $x = 0; $x <= "$#matrix"; $x++) {
	$max_length = "$#{$matrix[$x]}";
}
print "\n\nMatrix size = $nb_seq sequences ; $max_length positions\n";


### création des séquences consensus
open (ALI_OUT, ">$ARGV[1]");			# ouverture du fichier de sortie
for (my $i = 0; $i <= $#matrix; $i++) {	# pour chaque séquence
	$list_contigs{$seq_names[$i]}=$species_names[$i];		# on mémorise le couple nom du contig => nom de l'espèce
	$chevauchement{$species_names[$i]}=1;					# compteur à 0
	$difference{$species_names[$i]}=0;						# compteur à 0
	$identity{$species_names[$i]}=0;						# compteur à 0
	print "\n$seq_names[$i]";
	if ( !exists $seq_consensus{$species_names[$i]} ) {		# si le consensus de cette espèce n'existe pas encore
			$seq_consensus{$species_names[$i]}=$sequences{$seq_names[$i]};	# la séquence du contig actuelle devient le consensus
			print "\tfirst";
	}
	elsif ( exists $seq_consensus{$species_names[$i]} ) { 	# si le consensus existe déjà
		@tmp_consensus=split('', $seq_consensus{$species_names[$i]});		# on stock le consensus de manière temporaire
		$seq_consensus{$species_names[$i]}='';									# puis on vide le consensus original
		for (my $position = 0; $position < $max_length; $position++) {			# pour chaque position de la matrice
			if ($matrix[$i][$position] =~ m/[-\?X\*]/ && $tmp_consensus[$position] =~ m/[-\?X\*]/) {		# si position actuelle manquante et position consensus manquante
				$seq_consensus{$species_names[$i]}=$seq_consensus{$species_names[$i]}.$tmp_consensus[$position];	# on garde l'état du consensus
			}
			elsif ($matrix[$i][$position] eq $tmp_consensus[$position]) {			# si identité position actuelle/consensus
				$seq_consensus{$species_names[$i]}=$seq_consensus{$species_names[$i]}.$tmp_consensus[$position];	# on garde l'état du consensus
				++$chevauchement{$species_names[$i]};
				++$identity{$species_names[$i]};
			}
			elsif ($matrix[$i][$position] !~ m/[-\?X\*]/									# si actuelle est renseignée
			&& $tmp_consensus[$position] !~ m/[-\?X\*]/								# et consensus est renseignée
			&& $matrix[$i][$position] ne $tmp_consensus[$position]) {				# et désaccord
				$seq_consensus{$species_names[$i]}=$seq_consensus{$species_names[$i]}.$tmp_consensus[$position];	# on garde l'état du consensus (arbitraire, mais prk pas)
				++$chevauchement{$species_names[$i]};
				++$difference{$species_names[$i]};
			}
			elsif ($matrix[$i][$position] =~ m/[-\?X\*]/								# si actuelle est manquantes
			&& $tmp_consensus[$position] !~ m/[-\?X\*]/) {							# et consensus est renseignée
				$seq_consensus{$species_names[$i]}=$seq_consensus{$species_names[$i]}.$tmp_consensus[$position];	# on garde l'état du consensus
			}
			elsif ($matrix[$i][$position] !~ m/[-\?X\*]/									# si actuelle est renseignée
			&& $tmp_consensus[$position] =~ m/[-\?X\*]/) { 							# et consensus est manquante
				$seq_consensus{$species_names[$i]}=$seq_consensus{$species_names[$i]}.$matrix[$i][$position];	# on garde l'état de l'actuelle
			}
		}
#		if ($chevauchement{$species_names[$i]} > 25 && $difference{$species_names[$i]} > n) {		# si longueur chevauchement > 25 ET differences > n (pour poser un seuil sur le nbre de diff dans l'overlap)
		if ($chevauchement{$species_names[$i]} > 25) {												# si longueur chevauchement > 25 -> we consider it a 'significant' overlap
			$percentid = ($identity{$species_names[$i]}*100/$chevauchement{$species_names[$i]});		# calcul du id% pour le contig en court
			print "\twarning\tcoverage of $chevauchement{$species_names[$i]}\t(nb_diff = ",$difference{$species_names[$i]},")\t(id = ",$percentid,"%)";
			$alerte{$species_names[$i]}=1;															# alerte existe (et passe à 1) !!!!!!!
		}
		else {																								# sinon: chevauchement considered as 'acceptable'
			$percentid = ($identity{$species_names[$i]}*100/$chevauchement{$species_names[$i]});		# calcul du id% pour le contig en court
			print "\tFUSION\t(coverage of $chevauchement{$species_names[$i]}\t(nb_diff = ",$difference{$species_names[$i]},")\t(id = ",$percentid,"%)";
		}
	}
}


### ré-écriture de l'alignement
open (NOTIFY, ">$ARGV[2]");								# ouverture du fichier where we precise if overlap or not.
foreach my $key (keys(%seq_consensus)) {			# pour chaque espece (les clés de seq_consensus sont les espèces)
	@ctgs=();
	foreach my $ctgs (keys(%list_contigs)) {				# pour chaque contigs de l'alignement
		if ($list_contigs{$ctgs} eq $key) {					# si le contig appartient à l'espèce en court
			$ctgs =~ s/>//;									# mise en forme de nom de la séquence
			push (@ctgs, $ctgs);								# on le liste 
		}
	}
	if ($#ctgs+1 > 1 && ! exists $alerte{$key} ) {		# si il y a plusieurs contigs d'une espèce mais pas d'alerte chevauchement
		print ALI_OUT ">$key\|consensus_",join("_",@ctgs),"\n$seq_consensus{$key}\n";
		print NOTIFY "fusion\n"; 							# notifier que cette spS/famF est "sans-chevauchement" = on va fusionner et la seq consensus est dans output_name.out
	}
	elsif ($#ctgs+1 == 1) {								# ou si il y a un seul contig d'une espèce
		print ALI_OUT ">$ctgs[0]\n$seq_consensus{$key}";
		print NOTIFY "fusion\n"; 							# notifier que cette spS/famF est "sans-chevauchement" = on va fusionner et la seq consensus est dans output_name.out
	}
	elsif ($#ctgs+1 > 1 && exists $alerte{$key} ) {		# ou si il y a plusieurs contigs d'une espèce ET une alerte chevauchement (= overlap of more than 7 AA here)
		print NOTIFY "overlap\n"; 						# notifier que cette spS/famF est "avec-chevauchement" = on va devoir faire autrement
		#foreach my $ctgs2 (keys(%list_contigs)) {			# pour chaque contigs de l'alignement
			#if ($list_contigs{$ctgs2} eq $key) {				# si le contig appartient à l'espèce en court
				#my $ctgs3 = $ctgs2;
				#$ctgs2 =~ s/>//;										# mise en forme de nom de la séquence
				#print ALI_OUT ">$ctgs2\n$sequences{$ctgs3}";	# on la ré-écrit
				#print "\ntentative ecriture de :\n\tsp = $key\n\tctg = $ctgs3\n\tseq = $sequences{}";
			#}
		#}
		# optionnellement, on peut tout de même afficher le consensus des cas où l'on a pas le droit de fusionner
		#print ALI_OUT ">$key\@ALERTE-consensus_",join("_",@ctgs),"\n$seq_consensus{$key}\n";
	}
}

print "\n";
close ALI_OUT;
close NOTIFY;

# fin explicite
exit;

