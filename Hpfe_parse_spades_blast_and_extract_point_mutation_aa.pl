$file_base="D:\\antibiotic_resistance\\blastout/";
@file=("recA","rpsU","omp11","rdxA","FrxA","mdaB","sodB","fur","ribF");

foreach(@file){push @file_path_read, $file_base.$_.".blast";
			   push @file_write, $file_base.$_."_parse.txt";}

for(my $x=0;$x<=$#file_write;$x++){
	parse($file_path_read[$x],$file_write[$x]);
	rearrange($file_write[$x]);
	for_chi($file_write[$x]);
}

sub parse{
 my $file_path_read =shift @_;
 my $file_write=shift @_;
open FILE,"$file_path_read" or die "no such file $file_path_read";
local $/=">";
@blast=<FILE>;
close File;
shift @blast;
open File, ">$file_write";
for(my $a=0;$a<=$#blast;$a++){
	my @temp=split("\n",$blast[$a]);
	$temp[0]=~s/^\s+//;
	#print join("\n",@temp);
	push @head, $temp[0];
for(my $c=1;$c<=$#temp - 2;$c++){
	if($temp[$c]=~s/Query\s{2}(\d+)\s+//){
		$location=$1; #print $location,"\n";
		if($temp[$c+2]=~s/^Sbjct  \d+\s{2}//){		
		if(($temp[$c-3]=~/Identities = (\d+)\/\d+ \((\d+)%\)/) and $1 >60 and $2 > 90){
			#print "enter here \n";
			chomp($temp[$c]);
			chomp($temp[$c+1]);
			chomp($temp[$c+2]);
			$temp[$c]=~s/\s+(\d+)$//;
			$temp[$c+1]=~s/^\s+//;
			$temp[$c+2]=~s/\s+\d+$//;
			my @temp7=split(//,$temp[$c]);
			my @temp8=split(//,$temp[$c+1]);
			my @temp9=split(//,$temp[$c+2]);
			my $loc=$location-1;
			print File ">", $temp[0],"\n";
			for(my $b=0;$b<=$#temp7;$b++){
				$loc++;
				if($temp8[$b] eq " "  and $temp7[$b] ne $temp9[$b]){
					print File "mutation  $temp7[$b] $loc $temp9[$b]\n";
					if($temp7[$b] eq "-"){$loc--;}
				}
			}
		}
		$c=$c+2;
	 }
	}
}
}
close File;
}

## rearrange the mutations

sub rearrange{
	my $file_write=shift @_;
open File, $file_write;
my @file=<File>;
close File;
my @file2=split(/>/,join("",@file));
shift @file2;
#print @file2;
my %has_arr;
my %str_mut;
for(my $a=0;$a<=$#file2;$a++){
	my @temp=split "\n", $file2[$a];
	 $temp[0]=~/^\s*Hpfe(\d+)/;
	 $has_arr{$temp[0]}=$1;
	 my $count=$temp[0];
	 #print $count,"\n";
	 shift @temp;
	 $str_mut{$count}=join("\n", @temp);
}
open File, ">$file_write" or die "can not write";
foreach my $key (sort  { $has_arr{$a} <=> $has_arr{$b} } keys %has_arr){
	print File ">", $key,"\n",$str_mut{$key},"\n";
}
close File;
}
sub for_chi{
	my $file_write=shift @_;
	open File, $file_write;
	my @file=<File>;
	close File;
	my @file2=split(/>/,join("",@file));
	shift @file2;
	my %strain_mut;
	my %strain_hash_arr;
	my @mutation;
	my %temp;
	my %strain_info;
	
	### generate a pool of mutations
	for(my $a=0;$a<=$#file2;$a++){
		my @temp=split "\n", $file2[$a];
		$temp[0]=~/^(Hpfe\d+)/;
		$strain=$1;
		shift @temp;
		push @mutation, @temp;
		$strain_mut{$strain}=\@temp;
		$strain=~/Hpfe(\d+)/;
		$strain_hash_arr{$strain}=$1;
	}
	for(my $a=0;$a<=$#mutation;$a++){
		if($mutation[$a]=~/mutation  [-|A-Z|a-z] (\d+) [-|A-Z|a-z]/){
			$temp{$mutation[$a]}=$1;
		}else{ print "$file_write $strain $mutation[$a] not mutation pattern\n";}
	}
	@mutation =	sort  { $temp{$a} <=> $temp{$b} } keys %temp;
	
	#foreach $key(keys %strain_mut){ 
		#print $key, "\n",$#{$strain_mut{$key}},"\n"}
	#print $key,"\n", @{$strain_mut{$key}},"\n";
	#}
		
	#check the strains if mutations exists
	for(my $a=0;$a<=$#mutation;$a++){
		#print $mutation[$a],"list\n";
		foreach my $key(keys %strain_mut){
			#print $key,"\t",$strain_mut{$key},"\n";
			for(my $b=0;$b<=$#{$strain_mut{$key}};$b++){
				#print ${$strain_mut{$key}}[$b],"\n";
				if($strain_mut{$key}->[$b] eq $mutation[$a]){
					$strain_info{$key}->[$a]="y";
					last;
				}
				if(($b== $#{$strain_mut{$key}})  & ($strain_mut{$key}[$b] ne $mutation[$a])){
					$strain_info{$key}->[$a]="n";
				}
			}
		}
	}
	
	$file_write=$file_write.".tab";
	foreach(@mutation){s/mutation  //;}
	open File, ">$file_write" or die "can write file";
	print File "\t", join("\t", @mutation),"\n";
	foreach $key(sort {$strain_hash_arr{$a} <=> $strain_hash_arr{$b}} keys %strain_hash_arr){
		print File $key,"\t",join("\t", @{$strain_info{$key}}),"\n";
	}
	close File;
}
	