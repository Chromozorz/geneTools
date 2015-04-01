\p 5002
.g.ipc.h:0N;
.z.ps:{Q,::enlist x;if[.g.ipc.h;.g.ipc.h:hopen`::5001];loadvcf x};


loadvcf:{[path]
 vcfpath:hsym `$path;
 file:read0 hsym `$(path,"/","CEU.exon.2010_03.sites.vcf");
 metaend:first where not null {[file;x]if[(file x) like "#CHROM*";:x];}[file] each til count file;
 genome:("IJ*SSSS*";enlist "\t")0:(metaend)_file;
 formatinfo:4#'{ssr[;"-";""]each ssr[;"=";":`"]each x}each(";" vs '(exec INFO from genome));
 fmtcols:{[formatinfo;x]value"select ",("," sv formatinfo[x])," from (0#`)!()"}[formatinfo]each til count formatinfo;
 .g.ipc.genome:(delete INFO from genome),'fmtcols;
 neg[.g.ipc.h]("vcf sucessfully loaded. Table is called 'genome'");
 }
-1"This is the VCF engine";


