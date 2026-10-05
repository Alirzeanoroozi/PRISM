Bioinformatics, 2025, 41(12), btaf606
https://doi.org/10.1093/bioinformatics/btaf606
Advance Access Publication Date: 3 December 2025
Applications Note

Structural  bioinformatics

COCOMAPS 2.0: a web server for identifying, analyzing,
and visualizing atomic interactions at the interface of
biomolecular complexes
Mohit Chawla1,� , Utkarsh Kalra1, Andrea Petta2, Suraj Sharma3, Abdul Rajjak Shaikh4,
Luigi Cavallo1 , Romina Oliva4,�
1Physical Sciences and Engineering Division, King Abdullah University of Science and Technology (KAUST), Thuwal 23955-6900,
Saudi Arabia
2Tagetik Software S.R.L, Lucca 55100, Italy
3Department of Research and Innovation, STEMskills Research and Education Lab Private Limited, Faridabad 121002, Haryana, India
4Department of Sciences and Technologies, University “Parthenope” of Naples, Naples 80143, Italy
�Corresponding authors. Romina Oliva, Department of Sciences and Technologies, University “Parthenope” of Naples, Centro Direzionale Isola C4, Naples
80143, Italy. E-mail: romina.oliva@uniparthenope.it; Mohit Chawla, Physical Sciences and Engineering Division, King Abdullah University of Science and
Technology (KAUST), Thuwal 23955-6900, Saudi Arabia. E-mails: mohit.chawla@kaust.edu.sa; mohitchawla.bt@gmail.com.
Associate Editor: Jianlin Cheng

Abstract

Summary: Herein, we present COCOMAPS 2.0, for the analysis, visualization, and comparison of the interface in protein–protein and protein–
nucleic acid complexes. COCOMAPS 2.0 complements the residue-level and buried surface area analyses of the original COCOMAPS tool with
a  comprehensive  and  accurate  atomic-level  characterization  of  the  interface,  enabling  detailed  interpretation  of  molecular  recognition.
Furthermore, it provides a greatly enhanced flexibility, interactivity, and efficiency in graphical visualizations.

Availability and implementation: COCOMAPS 2.0 is accessible as a public web tool at https://aocdweb.com/BioTools/cocomaps2 and as a
standalone code at https://doi.org/10.5281/zenodo.17390665.

At  the  core  of  many  of  the  most  important  molecular  pro-
cesses  in  the  cell,  including  signal  transduction,  electron
transfer, gene expression and immune response, are interac-
tions between biomolecules. Increasing evidence is revealing
that  perturbation  of  such  interactions  frequently  leads  to
defective phenotypes. Mutations associated to human genetic
disorders tend to alter the interaction between proteins more
commonly than their folding and stability (Sahni et al. 2015,
Cheng  et  al.  2021,  Xiong  et  al.  2022).  Further,  disease-
related mutations in DNA-binding proteins have been shown
to  cluster  predominantly  at  the  DNA  interface,  where  they
affect  the  interaction  (Livesey  and  Marsh  2022,  Osterburg
et al. 2023).

An efficient visualization and analysis of the interface in the
3D structure of a biomolecular complex is therefore crucial for
a full  understanding of  the functional  and dysfunctional  bio-
logical processes driven by the associated interactions, as well
as for formulating testable predictions for interface modifica-
tion  and  targeting  (do  Nascimento  et  al.  2024,  Camps-Fajol
et al. 2025).

In  2011,  we  presented  COCOMAPS  (bioCOmplexes
COntact MAPS), a web tool using inter-residue contacts for
the  analysis  of  protein/nucleic  acid  complexes,  and  the  first
one  to  propose  inter-residue  contact  maps  for  an  efficient

visualization of the interface (Vangone et al. 2011). Also due
to its ease of use, COCOMAPS soon became popular and has
remained so. It has been used for instance for complementing
the  analysis  of  the  interface  in  newly  solved  experimental
structures of large complexes (Li et al. 2017, Whitehead et al.
2023, Zinzula et al. 2024), as well as for a large-scale compari-
son  of  SARS-CoV-2  spike  therapeutic  antibody  candidates
(Schendel et al. 2025), to cite recent applications.

In  the  last  decade,  availability  of  structures  of  assemblies
solved  experimentally  (Berman  et  al.  2000)  or  predictable
at  accuracy  competitive  with  experimental  structures  has
enormously increased (Lensink et al. 2023, Abramson et al.
2024).  This  motivated  us  to  provide  the  community  with  a
version 2.0 of COCOMAPS. In it, we put to use the under-
standing of the interface features we have been sharpening in
our multiple participations as scorers in the CAPRI (Critical
Assessment of PRedicted Interactions) experiment (Vangone
et al. 2013, Lensink et al. 2016, Lensink et al. 2017, Lensink
et  al.  2018,  Lensink  et  al.  2019,  Barradas-Bautista  et  al.
2020, Lensink et al. 2023), and combine it with our expertise
in  detecting  and  energetically  characterizing  non-covalent
interactions  in  biomolecules  (Chawla  et  al.  2014,  Chawla
et  al.  2015,  Chawla  et  al.  2017,  Kalra  et  al.  2020,  Chawla
et al. 2022).

Received: 2 June 2025; Revised: 20 October 2025; Accepted: 30 October 2025
© The Author(s) 2025. Published by Oxford University Press.
This is an Open Access article distributed under the terms of the Creative Commons Attribution License (https://creativecommons.org/licenses/by/4.0/), which
permits unrestricted reuse, distribution, and reproduction in any medium, provided the original work is properly cited.

Chawla et al.
2

Major advancements of COCOMAPS 2.0, as compared to
COCOMAPS,  include:  (i)  accepting  input  structures  in  the
mmCIF  file  format;  (ii)  detecting  atomic  interactions  and
accurately classifying them in 16 distinct classes; (iii) provid-
ing  an  interactive  Mol� visualization  (Sehnal  et  al.  2021),
which  encompasses  all  the  atomic  interactions  at  the  inter-
face, and features hyperlinks for seamless navigation between
3D view and tables; (iv) supplying multiple representations of
the interface, as a pie chart, a heatmap, a 2D and a 3D atomic
contact map.

The  core  of  the  COCOMAPS  2.0  code  lies  in  detecting
atomic  interactions  at  the  interface,  offering  what  is,  to  the
best of our knowledge, the most comprehensive classification
currently  available  among  tools  performing  similar  analyses
(Jubb  et  al.  2017,  Kayikci  et  al.  2018,  Badaczewska-Dawid
et al. 2022, Szulc et al. 2022, Del Conte et al. 2024, Schake
et  al.  2025).  A  comparative  summary  of  the  main  features
supported by such tools is provided in Table S1. The 16 types
of  interactions  we  classify  are:  H-bonds,  salt  bridges,  weak
H-bonds  (CH-ON  bonds),  water-mediated  contacts,  metal-
mediated contacts, disulfide bonds, halogen bonds, π–π/lone
pair-π/anion-π/cation-π/amino-π/ONSH-π/CH-π  interactions,
polar  and  apolar  vdW  contacts.  Most  of  these  interactions
are  variably  present  across  analogous  tools,  whereas  metal-
mediated, lone pair-π and amino-π contacts are not currently
defined in any of them (see Table 1, available as supplemen-
tary  data at  Bioinformatics  online).  Additionally,  contacts
with  interatomic  distances  below  the  sum  of  the  respective
vdW  radii  are  classified  as  ‘clashes’;  in  case  they  meet  the
criteria  for  one  of  the  interaction  types  listed  above,  this
interaction type is reported and supplemented with an asterisk
(�). Finally, contacts that fall within the selected cutoff inter-
residue distance (5 Å by default), but do not correspond to any
of  the  defined  atomic  contact  classes,  are  referred  to  as
‘proximal’ (Jubb et al. 2017).

Classification  of  the  interactions  is  based  on  criteria  we
carefully derived from the literature, including our own stud-
ies  (Chawla  et  al.  2014,  Chawla  et  al.  2017,  Kalra  et  al.
2020, Chawla et al. 2022). Programs under the COCOMAPS
2.0 web tool have been written in Python 3.0, taking advan-
tage of the numpy (v1.26.4), pandas (v2.2.2), scipy (v1.14.0),
biopython  (v1.83),  and  cctbx-base  (v2024.5)
libraries.
Reduce (Word et al. 1999) is used for adding hydrogen atoms
and HBPLUS (McDonald and Thornton 1994) for detecting
hydrogen  bonds,  salt  bridges,  and  water-mediated  contacts.
Other interactions are identified using in-house Python scripts.
The selected parameters and thresholds were validated in pre-
vious studies and are reported in detail, with the corresponding
references, on the About page of the web server (https://aocd
web.com/BioTools/cocomaps2/about).  Accessible  surfaces  are
calculated by NACCESS (Hubbard and Thornton 1993).

Inputs for COCOMAPS 2.0 are 3D structures, experimen-
tal or predicted, of complexes between protein, DNA and/or
RNA molecular chains. Input files, in the PDB or the mmCIF
format, can be directly retrieved from the wwPDB (with no
size limit; extended PDB IDs are also supported) or uploaded
locally  (size  up  to  40  MB).  Once  a  structure  has  been
uploaded,  the  input  page  of  COCOMAPS  2.0  will  list  the
chain IDs for the molecules present in it and their respective
range of residues and will allow the selection of the molecules
(or part of them) involved in the interaction to be analyzed.
In the advanced options, users may optionally also specify: a

project name, a name for the molecules involved in the inter-
action  and  their  email  address,  where  to  receive  the  output
once the job is done. Users can also modify the cutoff distance
for  defining  the  residue-residue  contacts  and  change,  within  a
given range, thresholds for all the geometrical parameters (dis-
tances and angles) used to detect and classify the atomic interac-
tions.  After  clicking  the  «Submit»  button,  users  are  redirected
to the output page, which checks the job status and reloads until
completion. Details on the runtime performance and scalability
are also given on the web server About page.

COCOMAPS 2.0 outputs are displayed on the results web
page for 3 months and archived as downloadable compressed
files. A link to the online resource is also emailed to the user,
if  requested.  At  the  top  of  the  COCOMAPS  2.0  output,  a
Table  of  Atomic  Interactions  and  a  Mol� 3D  view  of  the
complex,  with  interacting  residues  in  a  ball-and-stick  repre-
sentation, are shown side by side (see Fig. 1). Each row of the
Table  reports  a  pair  of  interacting  residues  and  lists  all
the  atomic  interactions  between  them.  The  table  is  sortable
by  residue  number  and  filterable  by  interaction  name/type.
A full interactivity for seamless navigation between the Table
and the 3D view is provided. When clicking a Table row, the
corresponding residues will be zoomed in and highlighted in
the  Mol� visualization.  By  clicking  the  «þ»  button  at  the
start of a row, all the atomic interactions between that pair
of  residues  will  be  listed  and  detailed  with  names  of  the
interacting  atoms  and  relative  distance/angle.  By  clicking
the eye-shaped button next to an interaction type, the corre-
sponding  interactions  will  be  zoomed  and  visualized  as
dashed lines in the Mol� visualization; the tag-shaped but-
ton activates labels on the interacting residues. On the other
hand,  by  clicking  on  a  residue  in  Mol�,  the  Table  will  be
sorted to display at the top and highlight all the interactions
involving that residue.

As  an  example,  in  Fig.  1,  the  COCOMAPS  2.0  3D  view
and Table of Atomic Interactions are shown for the barnase–
barstar complex (PDB ID: 1x1u), after clicking on the barstar
residue  Asp39  (chain  D)  at  the  interface.  The  table  appears
sorted to top list, highlighted in pale green, all the contacts in-
volving Asp39. These are 11 specific atomic interactions with
six different barnase residues, including 4 H-bonds, of which
3 are salt bridges, and 2 water-mediated contacts. For clarity,
only the 3 salt bridge interactions have been shown in the 3D
view and only that with Arg87 has been labeled. Both from
the table and for the Mol� view, it is immediately clear that
Asp39 is a hotspot for recognition. It is in fact the residue
contributing  most  to  the  complex  stability,  with  its  single
mutation  to  alanine  causing  a  drop  in  the  free  energy  of
binding  (ΔG)  by  7.7 kcal/mol  (Schreiber  and  Fersht  1995).
For  the  sake  of  comparison,  the  COCOMAPS  1.0  output
also reported the six Asp139 inter-residue contacts with bar-
nase;  however,  among  the  11  atomic  interactions,  only
the  4  hydrogen  bonds were  listed,  with no  specification  of
salt  bridges  (see  Fig.  1,  available  as  supplementary  data at
Bioinformatics online).

The Mol� 3D visualization presents, in the top right corner,
a selected Menu designed to resemble that in the «Explore in
3D»  view  of  the  wwPDB  (Berman  et  al.  2000),  to  simplify
the visualization experience to users already familiar with the
wwPDB. It includes a «Screenshot/State snapshot» button, for
taking high-definition pictures of the 3D view currently on the
screen, a «Toggle expanding viewpoint button», for full-screen

COCOMAPS 2.0

3

Figure 1. Details of the COCOMAPS 2.0 output for the interactions involving the barstar residue Asp39 in the barnase–barstar complex (PDB ID: 1x1u).
Left: Table of Atomic Interactions; for the sake of readability the column reporting the Asp39 residue is not displayed here. Right: Mol� view; the three
salt bridge interactions have been displayed and the interaction with Arg87 labeled. Table and Mol� star visualization are shown side by side, as they
appear in the server output page.

visualization of the complex, and Toggles for «Control panels»,
«Settings/Control  info»  and  «Selection  mode»,  from  which
experienced  users  can  access  all  the  advanced  Mol� func-
tionalities, such as sequence view and molecular component
selection,  modified  molecular  representation,  high-quality
rendering using lighting and outlines, plus additional geometric
measurements (distances, angles, dihedrals) and labeling. Above
the  Mol� view,  a  «Download»  button  allows  retrieval  of  the
analysis results, while a «Reset» button enables users to return
to the initial settings at any moment with a single click.

By scrolling down the output page, two tables are displayed,
summarizing  the  interface.  One  table  reports  the  number  of
residue-level interactions, while the other lists atomic interac-
tions  categorized  by  type.  A  hyperlink  to  the  Mol� 3D  view
allows users to visualize all interactions of a selected type di-
rectly within the molecular representation of the complex.

Above these Summary tables, four clickable symbols connect
to as many plots, providing a graphical overview of the inter-
face. These are a pie chart, a heatmap, a 2D and a 3D contact
map  of  the  atomic  interactions  (see  Fig.  2).  For  all  of  them,
ready-to-print, high-resolution images, in the PNG, JPEG, PDF
and SVG formats are provided. By hovering over the heatmap,
details  of  the  corresponding  interactions  are  visualized.  In  the
contact maps, all the specific types of interactions are displayed
with different symbols and can be filtered in/out; the third di-
mension of the 3D map is in fact represented by the interaction
type.  All  these  graphical  outputs  are  also  interactive  with  the
tabular  data  and  the  Mol� 3D  view.  At  the  bottom  of  the
Summary tables/plots, an ASA Table reports information about
the buried surface area upon complex formation, also detailed
per residue, for each involved molecule.

In Fig. 2, the pie-chart, heatmap and 3D contact map sum-
marizing the interface features are shown, again, for the bar-
nase–barstar  complex.  The  plots  highlight  the  variety  of
atomic interactions at the interface and the large number of
residues  involved.  It  can  also  be  noted  that  the  interface  is
rich in strong electrostatic interactions.  In fact, salt bridges,
H-bonds  and  water-mediated  contacts  represent  altogether
almost one third of the interactions, accounting for the tight

binding which makes this complex a model system for high-
affinity protein–protein interactions (Lee and Tidor 2001).

To  further  illustrate  the  features  and  scope  of  COCOMAPS
2.0,  we  also  applied  it  to  the  complex  between  Escherichia  coli
tRNACys  and  cysteinyl-tRNA  synthetase  (CysRS).  On  the
input  page,  for  the  tRNA  chain,  we  selected  only  the  anticodon
residues  (34–36)  and,  for  the  protein  chain,  the  C-terminal
anticodon-binding  domain  residues  (403–461).  The  anticodon  is
in fact expected to play a major role in recognition, as sub-
stitutions within it, particularly at G34, are known to dra-
matically affect the CysRS cysteinylation (Komatsoulis and
Abelson  1993).  The  Summary  of  interactions  table  and  a
Mol� 3D view from the COCOMAPS 2.0 output are shown
in Fig. 3. Because of the residue range selection made, inter-
actions  in  the  table  are  specific  to  the  anticodon.  They  in-
clude  15  directional  atomic  interactions  (7  H-bonds,  3
water-mediated  contacts,  3  CH–O/N  bonds,  1  π–π  and  1
CH–π contact), 8 of which involve G34. G34 interacts with
four CysRS residues, through 2 H-bonds, 3 water-mediated
contacts, 1 CH–O/N bond, 1 π–π  and 1 CH–π  contact (see
Table 2, available as supplementary data at Bioinformatics
online  for  details).  The  water-mediated  contacts,  shown  in
the 3D view of Fig. 3, are especially relevant as they create a
hydration  network  involving  both  Asp436  and  Arg423.
Previously,  only  the  G34  H-bonds  with  CysRS  Arg427  and
Asp436 and the π–π  contact with Trp432 had been described
(Hauenstein  et  al.  2004).  To  the  best  of  our  knowledge,  no
other public web server can both detect and interactively dis-
play in 3D all the atomic interactions occurring between G34
and  CysRS  (Table  2,  available  as  supplementary  data at
Bioinformatics online). This example thus helps illustrate how
COCOMAPS 2.0 expands the scope of similar tools by provid-
ing a more complete picture of molecular recognition, based on
both strong and weak interactions, and on the fundamental hy-
dration network, while requiring minimal effort from the user
to extract and visualize desired information.

In  conclusion,  as  shown  above,  COCOMAPS  2.0  is
designed  to  offer  an  intuitive,  informative  and  productive
experience, regardless of the user’s familiarity with biomo-
lecular  structures  or  the  chemical  nature  of  interactions.

Chawla et al.
4

Figure 2. COCOMAPS 2.0 graphical overviews of the interface for the barnase–barstar complex (PDB ID: 1x1u). Top: Pie chart of the atomic interactions.
Middle: Heatmap reporting the number of atomic interactions for each pair of interacting residues. Bottom: 3D atomic contact map: each type of interaction
is displayed with a different symbol on a separate layer; hovering over the symbols, the identity of residues involved in the interaction is displayed.

The accuracy and comprehensiveness of the atomic interac-
tions  classification,  combined  with  full  interactivity  be-
tween  tabular  information,  graphical  outputs  and  3D
visualization,  make  COCOMAPS  2.0  useful  for  a  wide

range  of  users—from  undergraduate  students  exploring
molecular interfaces for the first time to experienced struc-
tural biologists conducting detailed and informed analyses
of physico-chemical interactions.

COCOMAPS 2.0

5

Figure 3. COCOMAPS 2.0 output for the interaction between the tRNACys GCA anticodon and the CysRS C-ter domain (PDB ID: 1u0b). Top: Mol� view
of the complex, where only the water-mediated contacts, all involving G34, are displayed. Bottom: Summary of interactions table. The water-mediated
contacts can also be displayed in Mol� with a simple click on the eye-shaped button here.

Author contributions
Mohit  Chawla  (Conceptualization  [equal],  Data  curation
[equal],  Formal  analysis  [equal],  Investigation  [equal],
Methodology [equal], Project administration [equal], Software
[equal],  Supervision  [equal],  Validation [equal],  Visualization
[equal], Writing—original draft [equal]), Utkarsh Kalra (Data
curation  [equal],  Formal  analysis  [equal],  Investigation
[equal],  Software  [equal],  Validation  [equal],  Visualization
[equal]), Andrea Petta (Data curation [equal], Formal analy-
sis [equal], Methodology [equal], Software [equal], Validation
[equal],  Visualization  [equal]),  Suraj  Sharma  (Data  curation
[equal],  Formal  analysis  [equal],  Methodology  [equal],
Software [equal], Validation [equal]), Abdul Rajjak Shaikh
(Data  curation  [equal], Formal  analysis  [equal],  Investigation
[equal], Methodology [equal], Validation [equal], Visualization
[equal]),  Luigi  Cavallo  (Conceptualization  [equal],  Formal
analysis [equal], Funding acquisition [equal], Project adminis-
tration  [equal],  Supervision  [equal],  Visualization  [equal],

Writing—review  &  editing  [equal]),  and  Romina  Oliva
(Conceptualization  [equal],  Data  curation  [equal],  Formal
analysis  [equal],  Funding  acquisition  [equal],  Investigation
[equal],  Methodology  [equal],  Project  administration  [equal],
Supervision  [equal],  Validation [equal],  Visualization  [equal],
Writing—original  draft  [equal],  Writing—review  &  edit-
ing [equal]).

Supplementary data
Supplementary data are available at Bioinformatics online.

Conflict of interest: None declared.

Funding
This work was supported by the KAUST [URF/1/4384-01-01
and URF/1/4701-01-01 to L.C.] and Ministero dell’Universit�a

Chawla et al.
6

e della Ricerca, MUR [NextGeneration EU PRIN 2022 grant
2022HREZJT to R.O.].

Code availability
COCOMAPS 2.0 is accessible as a public web tool at https://
aocdweb.com/BioTools/cocomaps2 and as a standalone code
at https://doi.org/10.5281/zenodo.17390665.

References

Abramson J, Adler J, Dunger J et al. Accurate structure prediction of
biomolecular  interactions  with  AlphaFold  3.  Nature  2024;630:
493–500.

Badaczewska-Dawid AE, Nithin C, Wroblewski K et al. MAPIYA con-
tact  map  server  for  identification  and  visualization  of  molecular
interactions in proteins and biological complexes. Nucleic Acids Res
2022;50:W474–82.

Barradas-Bautista D, Cao Z, Cavallo L et al. The CASP13-CAPRI tar-
gets as case studies to illustrate a novel scoring pipeline integrating
CONSRANK  with  clustering  and
interface  analyses.  BMC
Bioinformatics 2020;21:262.

Berman  HM,  Westbrook  J,  Feng  Z  et  al.  The  protein  data  bank.

Nucleic Acids Res 2000;28:235–42.

Camps-Fajol C, Cavero D, Minguill�on J et al. Targeting protein–protein
interactions in drug discovery: modulators approved or in clinical
trials for cancer treatment. Pharmacol Res 2025;211:107544.

Chawla M, Abdel-Azeim S, Oliva R et al. Higher order structural effects
stabilizing the reverse Watson-Crick Guanine-Cytosine base pair in
functional RNAs. Nucleic Acids Res 2014;42:714–26.

Chawla M, Chermak E, Zhang Q et al. Occurrence and stability of lone
pair-pi  stacking  interactions  between  ribose  and  nucleobases  in
functional RNAs. Nucleic Acids Res 2017;45:11019–32.

Chawla M, Kalra K, Cao Z et al. Occurrence and stability of anion-pi
interactions between phosphate and nucleobases in functional RNA
molecules. Nucleic Acids Res 2022;50:11455–69.

Chawla  M,  Oliva  R,  Bujnicki  JM  et  al.  An  atlas  of  RNA  base  pairs
involving modified nucleobases with optimal geometries and accu-
rate energies. Nucleic Acids Res 2015;43:6714–29.

Cheng F, Zhao J, Wang Y et al. Comprehensive characterization of pro-
tein–protein interactions perturbed by disease mutations. Nat Genet
2021;53:342–53.

Del Conte A, Camagni GF, Clementel D et al. RING 4.0: faster residue
interaction  networks  with  novel  interaction  types  across  over
35,000  different  chemical  structures.  Nucleic  Acids  Res  2024;52:
W306–W312.

do Nascimento AM, Marques RB, Rold~ao AP et al. Exploring protein–
protein  interactions  for  the  development  of  new  analgesics.  Sci
Signal 2024;17:eadn4694.

Hauenstein S, Zhang C-M, Hou Y-M et al. Shape-selective RNA recog-
nition  by  cysteinyl-tRNA  synthetase.  Nat  Struct  Mol  Biol  2004;
11:1134–41.

Hubbard SJ, Thornton JM. 'NACCESS' Computer Program. London:
Department  of  Biochemistry  and  Molecular  Biology,  University
College, 1993.

Jubb  HC,  Higueruelo  AP,  Ochoa-Monta~no  B  et  al.  Arpeggio:  a  web
server for calculating and visualising interatomic interactions in pro-
tein structures. J Mol Biol 2017;429:365–71.

Kalra K, Gorle S, Cavallo L et al. Occurrence and stability of lone pair-
pi and OH-pi interactions between water and nucleobases in func-
tional RNAs. Nucleic Acids Res 2020;48:5825–38.

Kayikci M, Venkatakrishnan AJ, Scott-Brown J et al. Visualization and
analysis  of non-covalent contacts  using the protein  contacts  atlas.
Nat Struct Mol Biol 2018;25:185–94.

Komatsoulis GA, Abelson J. Recognition of tRNA(Cys) by Escherichia
coli cysteinyl-tRNA synthetase. Biochemistry 1993;32:7435–44.
Lee LP, Tidor B. Barstar is electrostatically optimized for tight binding

to barnase. Nat Struct Biol 2001;8:73–6.

Lensink MF, Brysbaert G, Nadzirin N et al. Blind prediction of homo-
and  hetero-protein  complexes:  the  CASP13-CAPRI  experiment.
Proteins 2019;87:1200–21.

Lensink  MF,  Brysbaert  G,  Raouraoua  N  et  al.  Impact  of  AlphaFold
on structure prediction of protein complexes: the CASP15-CAPRI
experiment. Proteins 2023;91:1658–83.

Lensink  MF,  Velankar  S,  Baek  M  et  al.  The  challenge  of  modeling
protein assemblies: the CASP12-CAPRI experiment. Proteins 2018;
86:257–73.

Lensink MF, Velankar S, Kryshtafovych A et al. Prediction of homoprotein
and  heteroprotein  complexes  by  protein  docking  and  template-based
modeling: a CASP-CAPRI experiment. Proteins 2016;84:323–48.

Lensink  MF,  Velankar  S,  Wodak  SJ.  Modeling  protein–protein  and
protein–peptide  complexes:  CAPRI  6th  edition.  Proteins  2017;
85:359–77.

Li X, Liu S, Jiang J et al. CryoEM structure of Saccharomyces cerevisiae
U1  snRNP  offers  insight  into  alternative  splicing.  Nat  Commun
2017;8:1035.

Livesey BJ,  Marsh  JA.  The  properties of  human  disease  mutations  at

protein interfaces. PLoS Comput Biol 2022;18:e1009858.

McDonald IK, Thornton JM. Satisfying hydrogen bonding potential in

proteins. J Mol Biol 1994;238:777–93.

Osterburg C, Ferniani M, Antonini D et al. Disease-related p63 DBD
mutations impair DNA binding by distinct mechanisms and varying
degree. Cell Death Dis 2023;14:274.

Sahni N, Yi S, Taipale M et al. Widespread macromolecular interaction
perturbations in human genetic disorders. Cell 2015;161:647–60.
Schake P, Bolz SN, Linnemann K et al. PLIP 2025: introducing protein–
protein interactions to the protein–ligand interaction profiler. Nucleic
Acids Res 2025;53:W463–5.

Schendel SL, Yu X, Halfmann PJ et al.; Coronavirus Immunotherapeutic
Consortium. A global collaboration for systematic analysis of broad-
ranging antibodies against the SARS-CoV-2 spike protein. Cell Rep
2025;44:115499.

Schreiber  G,  Fersht  AR.  Energetics  of  protein–protein  interactions:
analysis  of  the  barnase–barstar  interface  by  single  mutations  and
double mutant cycles. J Mol Biol 1995;248:478–86.

Sehnal D, Bittrich S, Deshpande M et al. Mol viewer: modern web app
for 3D visualization and analysis of large biomolecular structures.
Nucleic Acids Res 2021;49:W431–7.

Szulc NA, Mackiewicz Z, Bujnicki JM et al. fingeRNAt-a novel tool for
high-throughput analysis of nucleic acid-ligand interactions. PLoS
Comput Biol 2022;18:e1009783.

Vangone A, Cavallo L, Oliva R. Using a consensus approach based on
the  conservation  of  inter-residue  contacts  to  rank  CAPRI  models.
Proteins 2013;81:2210–20.

Vangone A, Spinelli R, Scarano V et al. COCOMAPS: a web application
to analyse and visualize contacts at the interface of biomolecular com-
plexes. Bioinformatics 2011;27:2915–6.

Whitehead  JD,  Decool  H,  Leyrat  C  et  al.  Structure  of  the  N-RNA/P
interface indicates mode of L/P recruitment to the nucleocapsid of
human metapneumovirus. Nat Commun 2023;14:7627.

Word JM, Lovell SC, Richardson JS et al. Asparagine and glutamine:
using hydrogen atom contacts in the choice of side-chain amide orien-
tation. J Mol Biol 1999;285:1735–47.

Xiong D, Lee D, Li L et al. Implications of disease-related mutations at
protein–protein interfaces. Curr Opin Struct Biol 2022;72:219–25.
Zinzula  L,  Beck  F,  Camasta  M  et  al.  Cryo-EM  structure  of  single-
layered  nucleoprotein–RNA  complex  from  Marburg  virus.  Nat
Commun 2024;15:10307.

© The Author(s) 2025. Published by Oxford University Press.
This is an Open Access article distributed under the terms of the Creative Commons Attribution License (https://creativecommons.org/licenses/by/4.0/), which permits
unrestricted reuse, distribution, and reproduction in any medium, provided the original work is properly cited.
Bioinformatics, 2025, 41, 1–6
https://doi.org/10.1093/bioinformatics/btaf606
Applications Note


