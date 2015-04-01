if[`.g in key `;delete from `.g];  
lf:{system"l /home/michael/data/geneTools.q"}
/#############################################Startup info################################################
welcome:{
 ("\n                                                                                                   ";
  "Hello and welcome to geneTools!! exit by typing 'exit'                                               ";
  "Useage can be found by entering 'genehelp'\n                                                         ";
  "NOTE: User input is restricted. To run sums etc then type 'revert'                                   ";
  "You can set the console to normal by starting the program with the -devmode command line option\n
 ")
 };

welcomeDev:{
 ("\n                                                                                                   ";
  "Hello and welcome to geneTools!!                                                                     ";
  "NOTE: You have chosen to run this program in devmode!\n
 ")
 };

genehelp:{
 ("\n                                                                                                   ";
  "See .g.funclist[] for a list of functions.                                                           ";
  "\nExample;\n\tTo see the code definition of the 'complement' function; '.g.funcs.complement[]'       ";
  "\nAs mentioned;aspects of the default q interface are disabled, only specified keywords can be used. ";
  "This can be changed by typing 'revert' which will reset the q interpreter to its defaults.           ";
  "This is active for a period of 30s. After 30s, it will be re-reverted so that the help features etc  ";
  "are again, enabled\n 
 ")
 };

.g.func.complement:{-1("\n                                                                               ";
               "This function returns the complementary strand of a positive strand DNA molecule.        ";
               "The output is 3-5. To run this call .g.comp[seq;`DNA] or .g.comp[seq;`RNA]\n
                      ");
 };

.g.funclist:`complement`gc`primers!("Computes the complement of DNA and RNA strings";
                                    "Computes the gc content of DNA and RNA strings";
                                    "Outputs primers given a DNA/RNA string and the region to be amplified");       
/#############################################Initialise##################################################
checkinputD:{
 x:-1_x;                              / take off line break char
 if[x like "genehelp";:-1 genehelp`]; / check for help
 if[x like "exit";exit 0];            / trying to exit?
 if[x like ".g*]";:.Q.s value x];     / allow exploration of .g
 if[x like "revert";                  / hand over to q  
    -1"Reversion active for 30 seconds";
    / reset input
    system"x .z.pi";
    / start time
    .geneutil.ftime:.z.T;
    timerevert[];
    system"t 100";
   ];
 }

timerevert:{[]
 .z.ts:{
  if[(.z.T- .geneutil.ftime)>00:00:30;
   .z.pi:checkinputD;
   system"t 0";
   -1"WARNING: Interpreter re-reverted! See 'genehelp' for details";
   ];
  };
 }
/#############################################Gene Utils##################################################
.g.util.paramchecker:{[params]
 if[`seq in key params;
    if[not 10h~type params`seq;:"Input sequence must be a string"];
   ];
 if[`mol in key params;
    if[not any(all each("DNA";"RNA")in(trim upper string params`mol));
       :"Input the type of molecule. e.g. `DNA. Your input was: `",raze string params`mol;
      ];
   ];
 }
/GC content calculation
.g.util.gc:{[seq;mol]
 .g.util.paramchecker[`seq`mol!(seq;mol)];
 $[10h~type seq;seq:`$'(seq);:"Input sequence must be a string"];
 tot:count seq;
 gc:count where seq in `G`C;
 gccont:(100%tot)*gc;
 :gccont;
 }

/basic annealing temp calculation
.g.util.ta:{[mol;seq]
 64.9+41*((.g.util.gc[seq;mol]%100)*(count seq)-16.4)%count seq
 }

/anhydrous molecular weight
.g.util.amw:{[seq;mol]
 if[not any(all each("DNA";"RNA")in(trim upper string mol));
    :"Input the type of molecule. e.g. `DNA. Your input was: `",raze string mol;
   ]; 
 amw:"(An x 313.21) + (Tn x 304.2) + (Cn x 289.18) + (Gn x 329.21) - 61.96";
 basecount:{[seq;base]count where (`$'seq)in base}[seq]each (`A`T`G`C;`A`U`G`C) (trim upper string mol)like "RNA"; 
 stop;
 }
/##########################################Complement calc################################################
.g.complement:{[seq;mol]
 dnamap:`A`T`C`G!`T`A`G`C;
 rnamap:`A`U`C`G!`U`A`G`C;
 $[10h~type seq;seq:`$'(seq);:"Input sequence must be a string"];
 if[(trim upper string mol)like "DNA";
    retseq:raze string dnamap seq;
   ]; 
 if[(trim upper string mol)like "RNA";
    retseq:raze string rnamap seq;
   ];
 :retseq;
 }
/##########################################Primer code ###################################################
.g.primers:{[seq;mol;startp;endp;dnac;saltc]
 ampseq:seq (startp-1)+til 1+endp-startp;
 sseq:$[startp<200;
        sseq:seq til startp;
        sseq:seq(endp-200)+til (endp-200)
       ];
 eseq:$[(endp+200)>count seq;
        endp:(endp)_seq;
        eseq:200#(endp)_seq
       ];
 :getprimers[mol;ampseq;sseq;eseq;dnac;saltc];
 }

getprimers:{[mol;queryseq;start;end;dnac;saltc]
 gccont:.g.util.gc[queryseq;mol];
 minl:15;
 maxl:30;
 lens:{[start;len]
  :raze{[start;minl;x]
    :(enlist(minl#(x _start)))!(enlist minl);
    }[start;len]each til (count start)-len;
  };
 allprimes:raze{
  flip select sequence:key x,length:value x from (0#`)!()
  }each lens[start]each minl+til(1+maxl-minl);
 :pickprimes[gccont;queryseq;allprimes;dnac;saltc;mol];
 }

pickprimes:{[gccont;queryseq;allprimes;dnac;saltc;mol]
 /salt adjusted calc
 adj:"Tm= 64.9 +41*(yG+zC-16.4)/(wA+xT+yG+zC)";
 taquery:.g.util.ta[mol;queryseq];
 /calculate annealing temperature for each primer candidate
 allprimes:update anneal:.g.util.ta[mol]each allprimes`sequence from allprimes;
 /filter primers so they don't have differenct annealing temps from the queryseq
 allprimes:select from allprimes where anneal within ((floor taquery)-1;(ceiling taquery)+1);
 if[0=count allprimes;:"No candidate primers found!"];
 }
/##########################################Oligo Calc code################################################
.g.oligocalc:{[seq;mol]
  
 }
/##########################################Fst code######################################################
.g.vcf:{[pathstring]
 .g.ipc.h:hopen`::5002;
 neg[.g.ipc.h](pathstring);
 }

/##########################################KickStartIt###################################################
$[`devmode in key .Q.opt .z.x;
 system"x .z.pi";
 .z.pi:checkinputD
 ];

$[`devmode in key .Q.opt .z.x;
  -1 welcomeDev`;
  -1 welcome`
 ];

\p 5001
.z.ps:{
 show x;
 if[(x((first x ss "genome")+til count "genome"))~"genome";
    genome::.g.ipc.h`.g.ipc.genome;
   ];
 system"x .z.ps"
 };
