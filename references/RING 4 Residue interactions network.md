Nucleic Acids Research , 2024, 52 , W306–W312
https://doi.org/10.1093/nar/gkae337
Advance access publication date: 30 April 2024
Web Server issue

RING 4.0: faster residue interaction networks with novel
int er action types across o v er 35,0 0 0 different chemical
structures
Alessio Del Conte  1 , Giorgia F. Camagni  1 , Damiano Clement el  1 , Gio vanni Minervini  1 ,
Alexander Miguel Monzon  2 , Carlo Fer rar i 2 , Damiano Piovesan  1 , * and Silvio C.E. Tosatto  1 , *
1 Department of Biomedical Sciences, University of Padova, Padova, Italy
2 Department of Information Engineering, University of Padova, Padova, Italy
* To whom correspondence should be addressed. Tel: +39 049 827 6269; Email: silvio.tosatto@unipd.it
Correspondence may also be addressed to Damiano Piovesan. Email: damiano.piovesan@unipd.it

Abstract
Residue interaction networks (RINs) are a valuable approach for representing contacts in protein str uct ures. RINs have been widely used in
various research areas, including the analysis of mutation effects, domain-domain communication, catalytic activity, and molecular dynamics
simulations. T he RING serv er is a po w erful tool to calculate non-co v alent molecular interactions based on geometrical parameters, providing
high-quality and reliable results. Here, we introduce RING 4.0, which includes significant enhancements for identifying both co v alent and non-
co v alent bonds in protein str uct ures. It now encompasses seven different interaction types, with the addition of π-h y drogen, halogen bonds and
metal ion coordination sites. The definitions of all a v ailable bond types ha v e also been refined and RING can now process the complete PDB
chemical component dictionary (o v er 350 0 0 different molecules) which pro vides atom names and co v alent connectivity inf ormation f or all kno wn
ligands. Optimization of the software has impro v ed e x ecution time b y an order of magnitude. The RING web server has been redesigned to
provide a more engaging and interactive user experience, incorporating new visualization tools. Users can now visualize all types of interactions
simultaneously in the str uct ure vie w er and netw ork component. T he w eb serv er, including e xtensiv e help and tutorials, is a v ailable from URL:
https:// ring.biocomputingup.it/ .

Gr aphical abstr act

Introduction

Proteins  are  large  macromolecules  whose  tridimensional
structure and motion in conjunction with the molecular en-
vironment, provide the functional toolkit of the cell. Non-
covalent interactions contribute to the stability of the folded
structures ( 1 ), they are responsible for the interaction with
other molecules in the formation of supra-molecular struc-
tures and are crucial in molecular communication such as al-
losterism and signaling processes, with the latter often medi-
ated by transient interaction events ( 2 ).

Network theory provides an alternative approach to de-
scribe protein structures emphasizing the non-covalent pro-
tein’s interactions. Residue interaction networks (RINs), also
known as protein structure networks (PSNs), represent amino
acids as nodes and the connections between them as edges
( 3 ). The definition of these connections, or contacts, typically
takes into account parameters such as physico-chemical prop-
erties and geometrical constraints (distances and angles), en-
abling the analysis of the network properties and characteri-
zation of specific types of interactions.

Received: March 8, 2024. Revised: April 9, 2024. Editorial Decision: April 16, 2024. Accepted: April 19, 2024
© The Author(s) 2024. Published by Oxford University Press on behalf of Nucleic Acids Research.
This is an Open Access article distributed under the terms of the Creative Commons Attribution-NonCommercial License
(https: // creativecommons.org / licenses / by-nc / 4.0 / ), which permits non-commercial re-use, distribution, and reproduction in any medium, provided the
original work is properly cited. For commercial re-use, please contact journals.permissions@oup.com

l

D
o
w
n
o
a
d
e
d

i

l

/

/

/

/

/

f
r
o
m
a
c
a
d
e
m
c
.
o
u
p
.
c
o
m
n
a
r
/
a
r
t
i
c
e
5
2
W
1
W
3
0
6
7
6
6
0
0
7
9
b
y
K
o
ç
U
n
v
e
r
s
i
t
y
S
u
n
a
K
r
a
ç
L
b
r
a
r
y
u
s
e
r
o
n
1
4
A
u
g
u
s
t

i

i

2
0
2
6

Nucleic Acids Research , 2024, Vol. 52, Web Server issue

W 307

Analysis of RINs is offered by a number of different solu-
tions ( 4 ), of particular success is the Residue Interaction Net-
work Generator (RING) which we have been developing and
maintaining since 2011 ( 5–7 ).

The analysis of RINs has progressed from initially demon-
strating  the  feasibility  of  network  parameters  in  static
structures  to  encompassing  a  range  of  more  specialized
applications.

In the last two years, RING has made significant contribu-
tions to numerous studies covering a diverse range of topics in
protein research, including binding and allosterism, interfaces,
protein stability, protein engineering, ligand binding predic-
tion, software development and molecular dynamics analysis
( 8–10 ).

A recent study examines the relationship between the size
of protein complex interfaces and the likelihood of cotrans-
lational assembly. RING was employed to calculate contact-
related properties of the interfaces, characterizing the different
kinds of interactions between subunits ( 11 ).

Structural ensemble analysis is one of the main applications
of RING. It is commonly used to analyze molecular dynam-
ics simulations by calculating probabilistic contact networks
( 12 ,13 ), and to analyze protein stability and perform virtual
screening ( 14 ,15 ).

Here, we present the new RING 4.0 version of the web
server that has been redesigned to offer a more engaging and
interactive user experience, including new visualization tools
and additional data. The underlying software library has been
optimized to enhance the calculation speed significantly. The
software has been improved and is now able to process all the
chemical compounds ( > 35 000 different molecules) available
in the PDB chemical component dictionary. RING is now able
to detect π-hydrogen and halogen bonds, as well as metal ion
coordination types in addition to previously available types,
whose definitions have been refined. RING 4.0 can now pro-
vide  comprehensive  interaction  networks  for  various  com-
binations of proteins, nucleotides, lipids, ligands, and small
molecules and all different types of interactions can be visu-
alized in the structure viewer and in a graph representation
simultaneously.

Materials and methods

The RING web server has at its core the RING application,
which has been extensively modified and improved compared
to the previous version ( 7 ). In terms of calculating interac-
tions, advancements were made to align with current stan-
dards. This involved broadening calculations to include any
ligand present in the PDB database, and introducing novel in-
teraction categories. Additionally, we improved the execution
time of about one order of magnitude. In the following, we
explain what changes have been introduced since the previ-
ous version, to the application and the web server.

Structure internal representation
The RING application takes as input protein structures, in the
format of PDB or PDBx / mmCIF file, requiring parsing and
conversion into an internal protein representation.

This internal representation aims to accurately reflect the
composition and chemical connectivity between atoms neces-
sary for interaction calculations. In the new version of RING
we expanded connectivity calculations beyond amino acids

and nucleic acids in order to detect typed interactions for all
possible compounds. Compounds chemical bonds and aro-
maticity are taken as external input from the PDB Chemical
Component Dictionary (CCD) ( 16 ) that contains over 35 000
entries and covers all structures available in the PDB. After es-
tablishing atom connections, a Depth First Search (DFS) algo-
rithm is employed to identify all possible rings, which are then
filtered based on the presence of all atoms listed as aromatic
in the CCD and their coplanarity within a fixed threshold.

RING can place hydrogen atoms in those input structures
lacking them, as the presence of hydrogens is crucial for many
interactions. When missing, RING can add hydrogens coor-
dinates based on the bond angle and dihedral angle with the
neighboring atoms. Since angles are fixed parameters and de-
pend on the type of atoms involved, in the new RING version
we retrieved all possible parameters for all compounds avail-
able in the PDB and stored them in a helper file.

During atom loading, various properties are assigned to
them. These properties include their van der Waals radius,
which has been expanded to encompass all known radii, as
well as indicators specifying whether an atom is a metal ion
or a halogen atom.

The integration of these data files into RING alleviates the
burden of real-time parameter calculations, mitigating the ex-
ecution time and possible errors. Moreover, users retain the
flexibility to modify this information as needed for customized
chemical components. When new compounds (not available
in the PDB) are analyzed, connectivity information, including
hydrogen placement, can be directly provided by the user.

Non-covalent bond calculation
In the latest version, significant modifications have been made
to  the  majority  of  prior  interactions,  including  hydrogen
bonding ( 17 ), π–π stacking ( 18 ), π–cation interactions ( 19 )
and van der Waals interactions. These changes were motivated
by several factors: expanding the definition of the interactions
to encompass all chemical components where feasible, opti-
mizing the software’s performance, and finally to align with
the latest advancements in defining non-covalent bonds as per
state-of-the-art publications.

Additionally, three new categories of interactions have been
introduced: the π–hydrogen bond ( 20 ), halogen bonds where
the acceptor can be an atom containing a lone pair of electrons
( 21 ) or a π system ( 22 ), and the metal ion coordination ( 23 ).
For a comprehensive understanding of how these interactions
are computed, please consult the documentation available on
the web server ( https:// ring.biocomputingup.it/ about ).

Execution time
In this updated version, significant improvements have been
made to the overall performance of the application. Previ-
ously, the implementation was designed to be versatile across
different uses, which meant the data structure was not well-
optimized for specific tasks. The latest update brings signifi-
cant enhancements to the software’s underlying logic and data
structures, resulting in a significant decrease in both memory
accesses and execution time. These improvements are facili-
tated by the implementation of efficient data structures like
hash maps, known for their average constant-time complex-
ity in operations such as search, insertion, and removal of el-
ements. This optimization expands the software’s utility be-
yond standard amino acids and nucleic acids, while also incor-

l

D
o
w
n
o
a
d
e
d

i

l

/

/

/

/

/

f
r
o
m
a
c
a
d
e
m
c
.
o
u
p
.
c
o
m
n
a
r
/
a
r
t
i
c
e
5
2
W
1
W
3
0
6
7
6
6
0
0
7
9
b
y
K
o
ç
U
n
v
e
r
s
i
t
y
S
u
n
a
K
r
a
ç
L
b
r
a
r
y
u
s
e
r
o
n
1
4
A
u
g
u
s
t

i

i

2
0
2
6

W 308

Nucleic Acids Research , 2024, Vol. 52, Web Server issue

porating the latest C++ standards to enrich the overall APIs.
Addressing performance issues, we utilized advanced profil-
ers like perf ( https://perf.wiki.kernel.org ), which allowed us
to identify and address time-consuming segments through it-
erative refinement of algorithms and optimization of memory
utilization, including the management of heap memory.

For instance, processing times for complex structures like
the RNA polymerase-Gfh1 (PDB code: 3aoh) have been no-
tably reduced, with RING 4 completing computations in just
3.2 s compared to the 22 s required by RING 3.0. Similarly,
handling larger structures such as the Virus-like particle of
bacteriophage AC (PDB code 6yf7) is now accomplished in a
mere 5.9 s, a considerable improvement over the 33 s previ-
ously needed. Even for conformational ensembles like the di-
hydrofolate reductase with multiple models (PDB code: 4p3r),
RING 4 demonstrates its efficiency, completing computations
in 9 s as opposed to the 30 s taken by RING 3.0.

Web-server implementation
The  web  server  utilizes  DRMAAtic  ( https://drmaatic.
biocomputingup.it )  as  its  backend,  an  in-house  solution
designed  for  executing  jobs  remotely  using  our  cluster  in-
frastructure.  DRMAAtic  interfaces  with  the  SLURM  ( 24 )
scheduler  via  the  Distributed  Resource  Manager  Applica-
tion  API  ( 25 ), managing  job  submission, permissions, user
management,  and  file  retrieval.  It  incorporates  throttling
mechanisms to prevent service misuse. Built on the Django
REST Framework ( 26 ) that provides most of the functionality
to create the REST APIs that are then used by the front-end
interface to submit the jobs and download the results. DR-
MAAtic facilitates the creation of REST APIs utilized by the
front-end  interface  for  job  submission  and  result  retrieval.
Authentication to DRMAAtic is achieved through OAuth 2.0
implicit flow, with ORCID ( https:// orcid.org/ ) serving as the
identity provider.

The RING w eb-serv er

Leveraging  on  DRMAAtic’s  capabilities,  the  web  server’s
front-end interface enables users to submit RING jobs with
customizable parameters and inputs, it will track the status of
the job, and retrieve the results when available. Users can login
to save submitted jobs to their account for easy retrieval across
different devices. The web server functions fully without re-
quiring user registration; past jobs are stored in the browser’s
cookies and can be associated with the account upon login.

Input
The home page of the web-server offers to the user an intro-
duction to what RING is, and an input form that can be used
for submitting new jobs. The form accepts as input a structure
file as PDB or PDBx / mmCIF formats, or an identifier that can
be both from the PDB ( 27 ) or AlphaFold DB ( 28 ). In the lat-
ter case, an autocomplete of the two has been implemented.
After selecting an input structure, the user can modify differ-
ent parameters, such as the chain to analyze (by default all
chains are considered), the model number (in case of multi-
ple state structure, all models are considered by default), and
many other options to change the interactions thresholds or
apply filters. While RING in principle has no limitations in the
number of models that can be processed, due to restrictions
in visualizing large structures, only files of up to 200 MB in

size are accepted by the server. For larger structural ensembles
we recommend installing RING locally . Finally , if the job sub-
mission passes all server checks, the user is redirected to the
result page, which checks the job status and reloads until it
completes. Otherwise an error message is displayed.

Output
When the RING job is finished, the front-end presents the out-
puts across various graphic components, to highlight different
aspects of the results. These include a 3D structure viewer built
with the Mol* Viewer ( 29 ), an interactive graph, representing
residues as nodes, and the interactions as edges and a con-
formational dependent contact map, available only when the
input has multiple states, shows for each residue if the selected
interaction is present or not in a state.

The new version of the output page has been improved in
several ways. First, the structure viewer has been enhanced to
reveal RING contacts and their associated details when fo-
cusing on specific residues on the graph. Second, we created a
probability contact map represented as a matrix highlighting
interacting residue pairs and their likelihood of contact across
multiple states.

Third, we integrated a contact table that lists interacting
residue pairs categorized by contact type, with probabilities
for multi-state inputs.

Notably, the structure viewer and the graph module can
be customized, it is possible to filter the interactions that are
displayed based on parameters like type, intra-chain or inter-
chain designation, and interaction frequency for multi-state
structures.

Users can also adjust residue color schemes based on prop-
erties such as chain, type, sequence position, secondary struc-
ture (DSSP) assignment, and node degree. Notably, all com-
ponents are tightly connected to each other enabling seamless
navigation; selecting a residue or interaction in one compo-
nent automatically highlights it across others, facilitating fo-
cused analysis for users.

The user can easily access and download various results,
including the original outputs generated by RING, its conver-
sion into a network saved in JSON file format, SVG or high-
resolution screenshots for other components, as well as tabu-
lar exports for the contact maps and tables.

Contact type distribution in PDB structures
We downloaded from the PDB (updated to February 2024)
all the structures with a resolution ≤1.6 Å, corresponding to
31 229 PDB structures. We executed RING on all these struc-
tures, and analyzed the distribution of contacts for different
classes of interactions. In Table 1 , we reported the distribu-
tion of contacts detected by RING, excluding results of PDB
entries without protein chains, therefore limiting the analy-
sis to 31 144 PDB structures, comprising 10 032 different
proteins (different UniProt identifiers). For proteins, we ex-
clusively consider standard amino acids, while for ligands, all
non-amino acid entities are included, and for nucleic acids,
both DNA and RNA bases are encompassed. Instead of using
the ‘auth_asym_id’ attribute, which corresponds to the chain
identifier in the canonical PDB format, we opted to consider
the value provided by the ‘label_asym_id’ attribute in the mm-
CIF format as chain identifiers. This decision was made as the
same chain identifier in the PDB format can be used to store
different molecular entities that are not necessarily covalently

l

D
o
w
n
o
a
d
e
d

i

l

/

/

/

/

/

f
r
o
m
a
c
a
d
e
m
c
.
o
u
p
.
c
o
m
n
a
r
/
a
r
t
i
c
e
5
2
W
1
W
3
0
6
7
6
6
0
0
7
9
b
y
K
o
ç
U
n
v
e
r
s
i
t
y
S
u
n
a
K
r
a
ç
L
b
r
a
r
y
u
s
e
r
o
n
1
4
A
u
g
u
s
t

i

i

2
0
2
6

Nucleic Acids Research , 2024, Vol. 52, Web Server issue

W 309

Table 1.  Distribution of RING contacts in high-resolution PDB str uct ures

Intra-chain

Inter-chain

Interaction partners

AA-AA

Number of PDBs
Number of contacts
Hydrogen bond
π- π stacking
π-cation
Ionic bond
Disulfide
π-hydrogen bond*
Halogen bond*
Metal ion coordination*

30 244
9 431 609

8 442 301
620 405
87 262
234 741
23 526
23 374
−
−

89.5%
6.6%
0.9%
2.5%
0.2%
0.2%
−
−

AA-Ligand ∧

3166
33 074

31 410
1475
102
−
−
35
43
−

95.0%
4.5%
0.3%
−
−
0.1%
0.1%
−

AA-AA

AA-Ligand

AA-Nucleic acid

12 605
300 687

229 295
36 301
7278
25 025
890
1898
−
−

76.3%
12.1%
2.4%
8.3%
0.3%
0.6%
−
−

25 062
324 296

229 855
33 175
2836
−
−
962
831
56 637

70.9%
10.2%
0.9%
−
−
0.3%
0.3%
17.5%

462
9311

8517
608
180
−
−
6
−
−

91.5%
6.5%
1.9%
−
−
0.1%
−
−

Only contacts involving at least one amino acid (AA) as an interacting partner are reported. Each column represents distinct interaction partners, and it is
indicated whether the interacting groups are in the same chain (Intra-chain) or different chains (Inter-chain). The ‘Number of PDBs’ refers to the total number
of PDB structures in which the interaction has been observed. For categorized interactions, we also provide the fraction of contacts observed for each category.
The statistics are derived by considering the ‘label_asym_id’ attribute in the mmCIF format as chain identifiers. (*) Interaction types introduced in the new
version of the RING software. ( ∧
) The intra-chain AA-Ligand column includes cases where standard amino acids interact with non-standard amino acids of
the same chain. Currently, non-standard amino acids are automatically processed by RING as ligands. However, in the current implementation of RING, ionic
bonds are detected only when amino acids are involved, and π-cation interactions are detected only when the cation is on the amino acid partner. Therefore,
their contribution in ligand interactions is missing.

bonded. By using the ‘label_asym_id’ attribute, we avoid the
potential pitfalls of inaccurate statistics and incorrect attribu-
tion of many inter-chain interactions to the intra-chain class.
Both the web server and the RING executable now offer the
option to choose which field to consider for the selection of
the chain identifier.

The analysis of non-covalent contacts across biomolecular

interactions in Table 1 unveils notable trends.

Hydrogen bonds emerge as ubiquitous as expected. The dif-
ference in the percentage of intra- and inter-chain AA–AA hy-
drogen bonds (89.5% and 76.3%) could be explained by the
different roles of this type of interaction in different contexts.
Intra-chain hydrogen bonds are predominant in the stabiliza-
tion of secondary structure elements, whereas surface inter-
chain interactions in protein-protein complexes are stabilized
by a more heterogeneous set of contacts that include a signif-
icant fraction of π–π stacking (12.1%) and π–cation (2.4%)
contacts among others.

π–π stacking  prominently  characterizes  both  protein–
protein (AA-AA) and protein–ligand (AA–ligand) inter-chain
interactions, indicating their importance in the stabilization of
supramolecular complexes.

Although  the  π–hydrogen  bond  is  relatively  uncommon
(0.3%, inter-chain AA–ligand), its presence in interactions be-
tween  proteins  and  ligands  emphasizes  its  potential  signif-
icance in stabilizing particular structural patterns. Protein–
ligand π–π stacking and π–hydrogen bonds are now detected
by the new version of RING for the first time thanks to the
fact that the software is now able to recognize aromatic rings
for all known PDB ligands.

As an example, RING was able to detect a very rare π–
hydrogen contact between an amino acid of the 50S ribosomal
protein L33 (PHE 39, label_asym_id / chain A) and a base of
a 23S rRNA (G 2374, label_asym_id / chain BA) of the 8b0a
PDB.

Metal ion coordination stands out prominently in ligand
interactions (AA–ligand) accounting for 17.5% of the total
inter-chain interactions. This highlights the significant role of
ionic coordination in molecular recognition and binding.

Finally, halogen bonds are extremely rare (0.3%) and only

observed in protein–ligand interactions.

Despite the great improvement in the classification of con-
tacts involving a ligand molecule, RING in its current imple-
mentation still presents some limitations. The software pro-
vides ionic bonds exclusively for amino acids, therefore ionic
contacts are not detected in protein–ligand complexes. More-
over, RING detects π–cation interactions only when the cation
is on the amino acid partner, therefore their abundance is un-
derestimated when a ligand is involved.

In conclusion, the three new types of interactions intro-
duced (metal ion coordination, π–hydrogen bond and halogen
bond) expand RING capacity of investigating processes such
as molecular recognition and protein–ligand binding.

Use case

The use case studies reported below illustrate how the new
RING  version  can  extract  valuable  biological  information
from molecular interactions.

The androgen receptor (AR) is a ligand-dependent tran-
scription  factor  involved  in  male  sexual  development  and
the etiology of pancreatic cancer (PCa). Therefore, it is con-
sidered  an  excellent  pharmacological  target  for  the  treat-
ment  of  prostate  cancer.  An  X-ray  crystallography  study
(PDB  code:  8e1a)  of  AR  in  complex  with  the  antagonist
compound VPC14368 is presented in Figure 1 A. The recep-
tor comprises three main functional domains: an N-terminal
domain  (NTD),  a  highly  conserved  DNA  binding  domain
(DBD),  and  a  C-terminal  ligand  binding  domain  (LBD).
VPC14368 binds to the androgen binding site (ABS) located in
the LBD.

RING generates a comprehensive network of interactions
between the protein and ligand, also expanding on intra-chain
protein-ligand interactions which were excluded in the pre-
vious version. This facilitates the study of molecular details
governing the AR-LBD ligand interaction (Figure 1 A). The
compound highlighted in the figure is stabilized in the bind-
ing pocket by a network of electrostatic interactions, includ-
ing van der Waals and aromatic interactions. The phenyl ring
of the ligand forms van der Waals interactions (gray) with
the surrounding side chains of M749, V746 and F764, estab-
lishing a π–π stack interaction (orange). The thiazole ring is

l

D
o
w
n
o
a
d
e
d

i

l

/

/

/

/

/

f
r
o
m
a
c
a
d
e
m
c
.
o
u
p
.
c
o
m
n
a
r
/
a
r
t
i
c
e
5
2
W
1
W
3
0
6
7
6
6
0
0
7
9
b
y
K
o
ç
U
n
v
e
r
s
i
t
y
S
u
n
a
K
r
a
ç
L
b
r
a
r
y
u
s
e
r
o
n
1
4
A
u
g
u
s
t

i

i

2
0
2
6

W 310

Nucleic Acids Research , 2024, Vol. 52, Web Server issue

l

D
o
w
n
o
a
d
e
d

i

l

/

/

/

/

/

f
r
o
m
a
c
a
d
e
m
c
.
o
u
p
.
c
o
m
n
a
r
/
a
r
t
i
c
e
5
2
W
1
W
3
0
6
7
6
6
0
0
7
9
b
y
K
o
ç
U
n
v
e
r
s
i
t
y
S
u
n
a
K
r
a
ç
L
b
r
a
r
y
u
s
e
r
o
n
1
4
A
u
g
u
s
t

i

i

2
0
2
6

Figure 1. Molecular insights into androgen receptor (AR) interactions calculated by RING 4.0. ( A ) AR-LBD complex with antagonist compound VPC14368
(PDB code: 8e1a) highlights electrostatic, van der Waals, and aromatic interactions. ( B ) AR’s DNA-binding domain (DBD) bound to DNA (PDB code: 1r4i)
sho w cases h y drogen bonding and π-cation interactions crucial f or DNA recognition.

stabilized by an electrostatic interaction with L704 and a π–π
stack with W741. Additionally, the morpholine group inter-
acts with the side chains of M780 and A877 via van der Waals
interactions (Figure 1 A).

As  a  transcription  factor,  the  androgen  receptor,  when
activated  by  endogenous  ligands  (testosterone  and  5 α-
dihydrotestosterone (DHT)), translocates into the nucleus and
initiates the transcription of target genes. In Figure 1 B, the
crystallized structure of the AR DNA-Binding Domain bound
to the Direct Repeat Response Element (PDB code: 1r4i) is
depicted. An interaction involving the R568 residue includes
two hydrogen bonds (blue) with guanine and a π-cation inter-
action (green) with the underlying selected nitrogenous base,
adenine. Furthermore, the α-helix, colored in purple, is stabi-
lized by π–π stack interactions between the nitrogenous bases.
RING 4.0 not only provides information on protein-protein
interactions  but  also  offers  molecular  details  on  protein-
ligand  and  protein-nucleic  acid  interactions. This  informa-
tion  can  be  utilized  to  study  differences  in  endogenous
protein-ligand and protein-inhibitor interactions or guide the
structure-based  design  of  new  molecules  targeting  specific
proteins.

more data. The help and about pages have been significantly
expanded to cover the new features.

This new version extends the previous server RING in sev-
eral ways. First of all, the underlying software library has been
carefully optimized to make the software about one order of
magnitude faster than RING 3.0 (or about 100 times faster
than RING 2.0).

Usage of the complete PDB chemical component dictionary
now provides the full range of 10 different non-covalent bond
types for over 35 000 different chemical compounds including
all known ligands. Of relevance, we have added π–hydrogen
and halogen bonds as well as metal ion coordination types.
Together with tightened definitions for all other non-covalent
bond types, RING 4.0 is now able to provide rich interac-
tion networks for any combination of protein, nucleotides,
lipids, ligands and small molecules. RING aims to deepen our
understanding of the intricate interplay of molecular forces
in  governing  biomolecular  stability  and  function,  offering
valuable insights for drug design and molecular engineering
endeavors.

The RING server is included in the Service Delivery Plan of
ELIXIR-Italy, implying a strong commitment to its continued
availability over the coming years.

Conclusions

RING has become a milestone for generating the increasingly
popular residue interaction networks. The RING web inter-
face has been re-designed to provide a richer and more inter-
active user experience, with novel visualization aids as well as

Data availability

No new data were generated or analysed in support of this re-
search. The web server, including extensive help and tutorials,
is available from URL: https:// ring.biocomputingup.it/ .

Nucleic Acids Research , 2024, Vol. 52, Web Server issue

W 311

Funding

ELIXIR,  the  research  infrastructure  for  life  science  data;
COST Action ML4NGP [CA21160], supported by COST (Eu-
ropean Cooperation in Science and Technology); European
Union’s Horizon 2020 research and innovation programme
[823886 MSCA-RISE ‘REFRACT’]; European Union through
NextGenerationEU,  PNRR  project  ELIXIRxNextGenIT
[IR0000010]; National Center for Gene Therapy and Drugs
based  on  RNA  Technology  [CN00000041];  Fondazione
AIRC  per  la  Ricerca  sul  Cancro  (AIRC)  [MFAG  29017
and IG 2019 ID 23825]; Italian Ministry of Education and
Research through the NextGenerationEU fund PRIN 2022
project:  PLANS  [2022W93FTW]. Funding  for  open  access
charge: University of Padova.

Conflict of interest statement

None declared.

References
1. Dill, K.A. and MacCallum, J.L. (2012) The protein-folding

problem, 50 years on. Science , 338 , 1042–1046.

2. del Sol, A. , Fujihashi, H. , Amoros, D. and Nussinov, R. (2006)
Residues crucial for maintaining short paths in network
communication mediate signaling in proteins. Mol. Syst. Biol., 2 ,
2006.0019.

3. del Sol,A. and Carbonell,P. (2007) The modular organization of
domain structures: insights into protein-protein binding. PLoS
Comput. Biol., 3 , e239.

4. Liang, Z. , Verkhivker, G.M. and Hu, G. (2020) Integration of

network models and evolutionary analysis into high-throughput
modeling of protein dynamics and allosteric regulation: theory,
tools and applications. Briefings Bioinf. , 21 , 815–835.

5. Martin, A.J.M. , V idotto, M. , Boscariol, F. , Domenico, T.D. , Walsh, I.

and Tosatto,S.C.E. (2011) RING: Networking interacting residues,
evolutionary information and energetics in protein structures.
Bioinformatics , 27 , 2003–2005.

6. Piovesan, D. , Minervini, G. and Tosatto, S.C.E. (2016) The RING
2.0 web server for high quality residue interaction networks.
Nucleic Acids Res., 44 , W367–W374.

7. Clementel, D. , Del Conte, A. , Monzon, A.M. , Camagni, G.F. ,

Minervini, G. , Piovesan, D. and Tosatto, S.C.E. (2022) RING 3.0:
fast generation of probabilistic residue interaction networks from
structural ensembles. Nucleic Acids Res. , 50 , W651–W656.
8. Verkhivker, G. , Agajanian, S. , Kassab, R. and Krishnan, K. (2022)

Probing mechanisms of binding and allostery in the S AR S-CoV-2
spike omicron variant complexes with the host receptor: revealing
functional roles of the binding hotspots in mediating epistatic
effects and communication with allosteric pockets. Int. J. Mol. Sci.,
23 , 11542.

9. Guo, H.-B. , Varaljay, V.A. , Kedziora, G. , Taylor, K. , Farajollahi, S. ,
Lombardo, N. , Harper, E. , Hung, C. , Gross, M. , Perminov, A. , et al.
(2023) Accurate prediction by AlphaFold2 for ligand binding in a
reductive dehalogenase and implications for PFAS (per- and
polyfluoroalkyl substance) biodegradation. Sci. Rep., 13 , 4082.
10. Krishnan, K. , T ian, H. , Tao, P. and Verkhivker, G.M. (2022) Probing

conformational landscapes and mechanisms of allosteric
communication in the functional states of the ABL kinase domain
using multiscale simulations and network-based mutational
profiling of allosteric residue potentials. J. Chem. Phys., 157 ,
245101.

11. Badonyi, M. and Marsh, J.A. (2022) Large protein complex

interfaces have evolved to promote cotranslational assembly.
eLife , 11 , e79602.

12. Sakai, T. , Mashima, T. , Kobayashi, N. , Ogata, H. , Duan, L. , Fujiki, R. ,
Hengphasatporn, K. , Uda, T. , Shigeta, Y. , Hifumi, E. , et al.  (2023)
Structural and thermodynamic insights into antibody light chain
tetramer formation through 3D domain swapping. Nat. Commun.,
14 , 7807.

13. Petrizzelli, F. , Biagini, T. , Bianco, S.D. , Liorni, N. , Napoli, A. ,

Castellana, S. and Mazza, T. (2022) Connecting the dots: a practical
evaluation of web-tools for describing protein dynamics as
networks. Front Bioinform , 2 , 1045368.

14. Ramakrishnan, K. , Johnson, R.L. , Winter, S.D. , Worthy, H.L. ,

Thomas, C. , Humer, D.C. , Spadiut, O. , Hindson, S.H. , Wells, S. ,
Barratt, A.H. , et al.  (2023) Glycosylation increases active site
rigidity leading to improved enzyme stability and turnover. FEBS
J., 290 , 3812–3827.

15. Kongsompong, S. , E-Kobon, T. , Taengphan, W. , Sangkhawasi, M. ,
Khongkow, M. and Chumnanpuen, P. (2023) Computer-aided
virtual screening and in vitro validation of biomimetic tyrosinase
inhibitory peptides from abalone peptidome. Int. J. Mol. Sci., 24 ,
3154.

16. Westbrook, J.D. , Shao, C. , Feng, Z. , Zhuravleva, M. , Velankar, S. and
Young,J. (2015) The chemical component dictionary: complete
descriptions of constituent molecules in experimentally determined
3D macromolecules in the Protein Data Bank. Bioinformatics , 31 ,
1274–1278.

17. Torshin, I.Y. , Weber, I.T. and Harrison, R.W. (2002) Geometric
criteria of hydrogen bonds in proteins and identification of
`bifurcated’ hydrogen bonds. Protein Eng. Des. Select., 15 ,
359–363.

18. Zhao, Y. , Li, J. , Gu, H. , Wei, D. , Xu, Y. , Fu, W. and Yu, Z. (2015)

Conformational preferences of π–π stacking between ligand and
protein, analysis derived from crystal structure data geometric
preference of π–π interaction. Interdiscip. Sci. Comput. Life Sci., 7 ,
211–220.

19. Infield, D.T. , Rasouli, A. , Galles, G.D. , Chipot, C. , Tajkhorshid, E. and
Ahern,C.A. (2021) Cation- π Interactions and their Functional
Roles in Membrane Proteins. J. Mol. Biol., 433 , 167035.
20. Steiner, T. and Koellner, G. (2001) Hydrogen bonds with

π-acceptors in proteins: frequencies and role in stabilizing local
3D structures. J. Mol. Biol., 305 , 535–557.

21. Chen, Z. , Wang, G. , Xu, Z. , Wang, J. , Yu, Y. , Cai, T. , Shao, Q. , Shi, J.
and Zhu,W. (2016) How do distance and solvent affect halogen
bonding involving negatively charged donors? J. Phys. Chem. B ,
120 , 8784–8793.

22. Shah, M.B. , Liu, J. , Zhang, Q. , Stout, C.D. and Halpert, J.R. (2017)
Halogen −π interactions in the cytochrome P450 active site:
structural insights into human CYP2B6 substrate selectivity. ACS
Chem. Biol., 12 , 1204–1210.

23. Zheng, H. , Cooper, D.R. , Porebski, P.J. , Shabalin, I.G. , Handing, K.B.

and Minor,W. (2017) CheckMyMetal: a macromolecular
metal-binding validation tool. Acta Crystallogr. D Struct. Biol., 73 ,
223–233.

24. Yoo, A.B. , Jette, M.A. and Grondona, M. (2003) SLURM: Simple

Linux Utility for Resource Management. In: Feitelson, D. ,
Rudolph, L. and Schwiegelshohn, U. (eds.) Job Scheduling Strategies
for Parallel Processing, Lecture Notes in Computer Science .
Springer Berlin Heidelberg, Berlin, Heidelberg, Vol. 2862 , pp.
44–60.

25. Troger, P. , Rajic, H. , Haas, A. and Domagalski, P. (2007)

Standardization of an API for distributed resource management
systems. In: Seventh IEEE International Symposium on Cluster
Computing and the Grid (CCGrid ’07) . pp. 619–626.

26. Tom Christie and Open source community (2022) Django REST

Framew or k . Version 3.14.

27. Burley, S.K. , Bhikadiya, C. , Bi, C. , Bittrich, S. , Chen, L. ,

Crichlow, G.V. , Christie, C.H. , Dalenberg, K. , Di Costanzo, L. ,
Duarte, J.M. , et al.  (2021) RCSB Protein Data Bank: powerful new
tools for exploring 3D structures of biological macromolecules for
basic and applied research and education in fundamental biology,

l

D
o
w
n
o
a
d
e
d

i

l

/

/

/

/

/

f
r
o
m
a
c
a
d
e
m
c
.
o
u
p
.
c
o
m
n
a
r
/
a
r
t
i
c
e
5
2
W
1
W
3
0
6
7
6
6
0
0
7
9
b
y
K
o
ç
U
n
v
e
r
s
i
t
y
S
u
n
a
K
r
a
ç
L
b
r
a
r
y
u
s
e
r
o
n
1
4
A
u
g
u
s
t

i

i

2
0
2
6

W 312

Nucleic Acids Research , 2024, Vol. 52, Web Server issue

biomedicine, biotechnology, bioengineering and energy sciences.
Nucleic Acids Res., 49 , D437–D451.

28. Varadi, M. , Anyango, S. , Deshpande, M. , Nair, S. , Natassia, C. ,
Yordanova, G. , Yuan, D. , Stroe, O. , Wood, G. , Laydon, A. , et al.
(2022) AlphaFold Protein Structure Database: massively
expanding the structural coverage of protein-sequence space with
high-accuracy models. Nucleic Acids Res. , 50 , D439–D444.

29. Sehnal, D. , Bittrich, S. , Deshpande, M. , Svobodová, R. , Berka, K. ,

Bazgier, V. , Velankar, S. , Burley, S.K. , Ko ˇca, J. and Rose, A.S. (2021)
Mol* Viewer: modern web app for 3D visualization and analysis
of large biomolecular structures. Nucleic Acids Res., 49 ,
W431–W437.

l

D
o
w
n
o
a
d
e
d

i

l

/

/

/

/

/

f
r
o
m
a
c
a
d
e
m
c
.
o
u
p
.
c
o
m
n
a
r
/
a
r
t
i
c
e
5
2
W
1
W
3
0
6
7
6
6
0
0
7
9
b
y
K
o
ç
U
n
v
e
r
s
i
t
y
S
u
n
a
K
r
a
ç
L
b
r
a
r
y
u
s
e
r
o
n
1
4
A
u
g
u
s
t

i

i

2
0
2
6

Received: March 8, 2024. Revised: April 9, 2024. Editorial Decision: April 16, 2024. Accepted: April 19, 2024
© The Author(s) 2024. Published by Oxford University Press on behalf of Nucleic Acids Research.
This is an Open Access article distributed under the terms of the Creative Commons Attribution-NonCommercial License
(https: // creativecommons.org / licenses / by-nc / 4.0 / ), which permits non-commercial re-use, distribution, and reproduction in any medium, provided the original work is properly cited. For
commercial re-use, please contact journals.permissions@oup.com


