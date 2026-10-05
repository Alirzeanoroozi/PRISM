proteins
STRUCTURE OFUNCTION OBIOINFORMATICS
Prediction of homoprotein and heteroprotein
complexes by protein docking and template-
based modeling: A CASP-CAPRI experiment
Marc F. Lensink,1* Sameer Velankar,2 Andriy Kryshtafovych,3 Shen-You Huang,4
Dina Schneidman-Duhovny,5,6 Andrej Sali,5,6,7 Joan Segura,8 Narcis Fernandez-Fuentes,9
Shruthi Viswanath,10,11 Ron Elber,11,12 Sergei Grudinin,13,14 Petr Popov,13,14,15
Emilie Neveu,13,14 Hasup Lee,16 Minkyung Baek,16 Sangwoo Park,16 Lim Heo,16 Gyu Rie Lee,16
Chaok Seok,16 Sanbo Qin,17 Huan-Xiang Zhou,17 David W. Ritchie,18 Bernard Maigret,19
Marie-Dominique Devignes,19 Anisah Ghoorah,20 Mieczyslaw Torchala,21
Rapha€el A.G. Chaleil,21 Paul A. Bates,21 Efrat Ben-Zeev,22 Miriam Eisenstein,23
Surendra S. Negi,24 Zhiping Weng,25 Thom Vreven,25 Brian G. Pierce,25 Tyler M. Borrman,25
Jinchao Yu,26 Franc¸oise Ochsenbein,26 Rapha€el Guerois,26 Anna Vangone,27
Jo~ao P.G.L.M. Rodrigues,27 Gydo van Zundert,27 Mehdi Nellen,27 Li Xue,27 Ezgi Karaca,27
Adrien S.J. Melquiond,27 Koen Visscher,27 Panagiotis L. Kastritis,27 Alexandre M.J.J. Bonvin,27
Xianjin Xu,28 Liming Qiu,28 Chengfei Yan,28,29 Jilong Li,30 Zhiwei Ma,28,29 Jianlin Cheng,30,31
Xiaoqin Zou,28,29,31,32 Yang Shen,33 Lenna X. Peterson,34 Hyung-Rae Kim,34 Amit Roy,34,35
Xusi Han,34 Juan Esquivel-Rodriguez,36 Daisuke Kihara,34,36 Xiaofeng Yu,37 Neil J. Bruce,37
Jonathan C. Fuller,37 Rebecca C. Wade,37,38,39 Ivan Anishchenko,40 Petras J. Kundrotas,40
Ilya A. Vakser,40,41 Kenichiro Imai,42 Kazunori Yamada,42 Toshiyuki Oda,42
Tsukasa Nakamura,43 Kentaro Tomii,42,43 Chiara Pallara,44 Miguel Romero-Durana,44
Brian Jim(cid:2)enez-Garc(cid:2)ıa,44 Iain H. Moal,44 Juan F(cid:2)ernandez-Recio,44 Jong Young Joung,45
Jong Yun Kim,45 Keehyoung Joo,45,46 Jooyoung Lee,45,47 Dima Kozakov,48 Sandor Vajda,48,49
Scott Mottarella,48 David R. Hall,48 Dmitri Beglov,48 Artem Mamonov,48 Bing Xia,48
Tanggis Bohnuud,48 Carlos A. Del Carpio,50,51 Eichiro Ichiishi,52 Nicholas Marze,53
Daisuke Kuroda,53 Shourya S. Roy Burman,53 Jeffrey J. Gray,53,54 Edrisse Chermak,55
Luigi Cavallo,55 Romina Oliva,56 Andrey Tovchigrechko,57 and Shoshana J. Wodak58,59*
AdditionalSupportingInformationmaybefoundintheonlineversionofthisarticle. VIB Structural Biology Research Center, VUB, 1050 Brussels, Belgium. E-mail:
Grant sponsor: NIH; Grant numbers: R01 GM083960; P41 GM109824; shoshana.wodak@gmail.com
GM058187; R01 GM061867; R01 GM093147; R01 GM078221; R01GM109980; Tyler M. Borrman current address is Institute for Bioscience and Biotechnology
R01GM094123; R01 GM097528; R01GM074255; Grant sponsor: Biotechnology Research,UniversityofMaryland,Rockville,MD20850.
and Biological Sciences Research Council; Grant number: BBS/E/W/10962A01D; Yang Shen current address is Center for Bioinformatics and Genomic Systems
Grantsponsor:ResearchCouncilsUKAcademicFellowshipprogram;Grantspon- Engineering, Department of Electrical and Computer Engineering, Texas a&M
sor:CancerResearchUK;Grantsponsor:KlausTschiraFoundation;Grantspon- University,CollegeStation,TX77843.
sor:PlatformProjectforSupportinginDrugDiscoveryandLifeScienceResearch; KenichiroImai,ToshiyukiOda,andKentaroTomiicurrentaddressisBiotechnol-
Grantsponsor:JapanAgencyforMedicalResearchandDevelopment;Grantspon- ogyResearchInstituteforDrugDiscovery,NationalInstituteofAdvancedIndus-
sor: Agence Nationale de la Recherche; Grant number: ANR-11-MONU-0006; trialScienceandTechnology(AIST),Koto-Ku,Japan.
Grant sponsor: National Research Foundation of Korea (NRF); Grant numbers: KazunoriYamadacurrentaddressisGroupofElectricalEngineering,Communica-
NRF-2013R1A2A1A09012229;2008-0061987;Grantsponsor:BIP;Grantnumber: tion Engineering, Electronic Engineering, and Information Engineering, Tohoku
ANR-IAB-2011-16-BIP:BIP;Grantsponsor:H2020MarieSklodowska-CurieIndi- University,Sendai,Japan.
vidualFellowship;Grantnumber:659025-BAP;Grantsponsor:NetherlandsOrga- IainH.MoalcurrentaddressisEuropeanMolecularBiologyLaboratory,European
nizationforScientificResearchVeni;Grantnumber:722.014.005;Grantsponsor: BioinformaticsInstitute(EMBL-EBI),WellcomeTrustGenomeCampus, Hinxton,
NationalScienceFoundation;Grantnumbers:CAREERAwardDBI0953839;CCF- CambridgeCB101SD,UnitedKingdom.
1546278;NSFIIS1319551;NSFDBI1262189;NSFIOS1127027;NSFDBI1262621; Daisuke Kuroda current address is School of Pharmacy, Showa University,
NSF DBI 1458509; NSF AF 1527292; Grant sponsor: EU; Grant number: FP7 Shinagawa-Ku,Tokyo142-8555,Japan.
604102 (HBP); Grant sponsor: BMBF; Grant number: 0315749 (VLN); Grant Thecopyrightlineforthisarticlewaschangedon6October2016afteroriginal
sponsor: Spanish Ministry of Economy and Competitiveness; Grant number: onlinepublication.
BIO2013-48213-R; Grant sponsor: European Union; Grant number: FP7/2007- ThisisanopenaccessarticleunderthetermsoftheCreativeCommonsAttribu-
2013REAPIEF-GA-2012-327899;Grantsponsor:NationalInstituteofSupercom- tion License, which permits use, distribution and reproduction in any medium,
puting and Networking; Grant number: KSC-2014-C3-01; Grant sponsor: US- providedtheoriginalworkisproperlycited.
Israel BSF; Grant number: 2009418; Grant sponsor: Regione Campania; Grant Received29May2015;Revised30December2015;Accepted2February2016
number:LR5-AF2008. Publishedonline28April2016inWileyOnlineLibrary(wileyonlinelibrary.com).
*Correspondence to: Marc F. Lensink; University Lille, CNRS UMR8576 UGSF, DOI:10.1002/prot.25007
Lille,F-59000,France.E-mail:marc.lensink@univ-lille1.frorShoshanaJ.Wodak;
VVC 2016THEAUTHORSPROTEINS:STRUCTURE,FUNCTION,ANDBIOINFORMATICSPUBLISHEDBYWILEYPERIODICALS,INC. PROTEINS 323

M.F.Lensinketal.
1UniversityLille,CNRSUMR8576UGSF,Lille,F-59000,France
2
EuropeanMolecularBiologyLaboratory,EuropeanBioinformaticsInstitute(EMBL-EBI),WellcomeTrustGenomeCampus,Hinxton,Cambridge,CB101SD,United
Kingdom
3
GenomeCenter,UniversityofCalifornia,Davis,California,95616
4
ResearchSupportComputing,UniversityofMissouriBioinformaticsConsortium,andDepartmentofComputerScience,UniversityofMissouri,
Columbia,Missouri65211
5
DepartmentofBioengineeringandTherapeuticSciences,UniversityofCaliforniaSanFrancisco,SanFrancisco,California94158
6
DepartmentofPharmaceuticalChemistry,UniversityofCaliforniaSanFrancisco,SanFrancisco,California94158
7
CaliforniaInstituteforQuantitativeBiosciences(QB3),UniversityofCaliforniaSanFrancisco,SanFrancisco,California94158
8
GN7oftheNationalInstituteforBioinformatics(INB)andBiocomputingUnit,NationalCenterofBiotechnology(CSIC),Madrid,28049,Spain
9
InstituteofBiological,EnvironmentalandRuralSciences(IBERS),AberystwythUniversity,Aberystwyth,SY233FG,UnitedKingdom
10
DepartmentofComputerScience,UniversityofTexasatAustin,Austin,Texas78712
11InstituteforComputationalEngineeringandSciences,UniversityofTexasatAustin,Austin,Texas78712
12DepartmentofChemistry,UniversityofTexasatAustin,Austin,Texas78712
13LJK,UniversityGrenobleAlpes,CNRS,Grenoble,38000,France
14
INRIA,Grenoble,38000,France
15
MoscowInstituteofPhysicsandTechnology,Dolgoprudniy,Russia
16
DepartmentofChemistry,SeoulNationalUniversity,Seoul,151-747,RepublicofKorea
17
DepartmentofPhysicsandInstituteofMolecularBiophysics,FloridaStateUniversity,Tallahassee,Florida32306,USA
18 INRIANancy—GrandEst,Villers-le`s-Nancy,54600,France
19 CNRS,LORIA,CampusScientifique,BP239,Vandœuvre-le`s-Nancy,54506,France
20
DepartmentofComputerScienceandEngineering,UniversityofMauritius,Reduit,Mauritius
21
BiomolecularModellingLaboratory,theFrancisCrickInstitute,Lincoln’sInnFieldsLaboratory,London,WC2A3LY,UnitedKingdom
22
G-INCPM,WeizmannInstituteofScience,Rehovot,7610001,Israel
23
DepartmentofChemicalResearchSupport,WeizmannInstituteofScience,Rehovot,7610001,Israel
24
SealyCenterforStructuralBiologyandMolecularBiophysics,UniversityofTexasMedicalBranch,301UniversityBoulevard,Galveston,Texas77555-0857
25
PrograminBioinformaticsandIntegrativeBiology,UniversityofMassachusettsMedicalSchool,Worcester,Massachusetts01605
26
InstituteforIntegrativeBiologyoftheCell(I2BC),CEA,CNRS,UniversityParis-Saclay,CEA-Saclay,Gif-sur-Yvette,91191,France
27
BijvoetCenterforBiomolecularResearch,FacultyofScience–Chemistry,UtrechtUniversity,Padualaan8,Utrecht,3584CH,TheNetherlands
28
DaltonCardiovascularResearchCenter,UniversityofMissouri,Columbia,Missouri65211
29
DepartmentofPhysicsandAstronomy,UniversityofMissouri,Columbia,Missouri65211
30
DepartmentofComputerScience,UniversityofMissouri,Columbia,Missouri65211
31
InformaticsInstitute,UniversityofMissouri,Columbia,Missouri65211
32
DepartmentofBiochemistry,UniversityofMissouri,Columbia,Missouri65211
33
ToyotaTechnologicalInstituteatChicago,6045SKenwoodAvenue,Chicago,Illinois60637
34
DepartmentofBiologicalSciences,PurdueUniversity,WestLafayette,Indiana47907
35
BioinformaticsandComputationalBiosciencesBranch,RockyMountainLaboratories,NationalInstitutesofHealth,Hamilton,Montano59840
36
DepartmentofComputerScience,PurdueUniversity,WestLafayette,IN,USA47907
37
MolecularandCellularModelingGroup,HeidelbergInstituteforTheoreticalStudies(HITS),Heidelberg,Germany
38
CenterforMolecularBiology(ZMBH),DKFZ-ZMBHAlliance,HeidelbergUniversity,Heidelberg,Germany
39
InterdisciplinaryCenterforScientificComputing(IWR),HeidelbergUniversity,Heidelberg,Germany
40
CenterforComputationalBiology,TheUniversityofKansas,Lawrence,Kansas66047
41
DepartmentofMolecularBiosciences,TheUniversityofKansas,Lawrence,Kansas66047
42
ComputationalBiologyResearchCenter(CBRC),NationalInstituteofAdvancedIndustrialScienceandTechnology(AIST),Koto-Ku,Japan
43
GraduateSchoolofFrontierSciences,theUniversityofTokyo,Kashiwa,Japan
44JointBSC-CRG-IRBResearchPrograminComputationalBiology,BarcelonaSupercomputingCenter,C/JordiGirona29,Barcelona,08034,Spain
45Centerforin-SilicoProteinScience,KoreaInstituteforAdvancedStudy,Seoul,130-722,Korea
46CenterforAdvancedComputation,KoreaInstituteforAdvancedStudy,Seoul,130-722,Korea
47SchoolofComputationalScience,KoreaInstituteforAdvancedStudy,Seoul,130-722,Korea
48
DepartmentofBiomedicalEngineering,BostonUniversity,Boston,Massachusetts
49
DepartmentofChemistry,BostonUniversity,Boston,Massachusetts
50
InstituteofBiologicalDiversity,InternationalPacificInstituteofIndiana,Bloomington,Indiana47401
51
DrosophilaGeneticResourceCenter,KyotoInstituteofTechnology,Ukyo-Ku,616-8354,Japan
52
InternationalUniversityofHealthandWelfareHospital(IUHWHospital),Asushiobara-City,TochigiPrefecture329-2763,Japan
53
DepartmentofChemicalandBiomolecularEngineering,JohnsHopkinsUniversity,Baltimore,Maryland21218
54
PrograminMolecularBiophysics,JohnsHopkinsUniversity,Baltimore,Maryland21218
55
KingAbdullahUniversityofScienceandTechnology,SaudiArabia
56
UniversityofNaples“Parthenope”,Napoli,Italy
57
J.CraigVenterInstitute,9704MedicalCenterDrive,Rockville,Maryland20850
58
DepartmentsofBiochemistryandMolecularGenetics,UniversityofToronto,Toronto,Ontario,Canada
59
VIBStructuralBiologyResearchCenter,VUBPleinlaan2,Brussels,1050,Belgium
324 PROTEINS

PredictionofHomoandHeteroproteinComplexesbyProteinDockingandModeling
ABSTRACT
We present the results for CAPRI Round 30, the first joint CASP-CAPRI experiment, which brought together experts from
the protein structure prediction and protein–protein docking communities. The Round comprised 25 targets from amongst
those submitted for the CASP11 prediction experiment of 2014. The targets included mostly homodimers, a few homo-
tetramers,andtwoheterodimers,andcomprisedproteinchainsthatcouldreadilybemodeledusingtemplatesfromthePro-
tein Data Bank. On average 24 CAPRI groups and 7 CASP groups submitted docking predictions for each target, and 12
CAPRI groups per target participated in the CAPRI scoring experiment. In total more than 9500 models were assessed
against the 3D structures of the corresponding target complexes. Results show that the prediction of homodimer assemblies
by homology modeling techniques and docking calculations is quite successful for targets featuring large enough subunit
interfaces to represent stable associations. Targets with ambiguous or inaccurate oligomeric state assignments, often featur-
ing crystal contact-sized interfaces, represented a confounding factor. For those, a much poorer prediction performance was
achieved, while nonetheless often providing helpful clues on the correct oligomeric state of the protein. The prediction per-
formance was very poor for genuine tetrameric targets, where the inaccuracy of the homology-built subunit models and the
smaller pair-wise interfaces severely limited the ability to derive the correct assembly mode. Our analysis also shows that
docking procedures tend to perform better than standard homology modeling techniques and that highly accurate models
of the protein components are notalwaysrequiredtoidentify their association modeswith acceptableaccuracy.
Proteins2016;84(Suppl1):323–348.
VC 2016TheAuthorsProteins:Structure,Function,andBioinformaticsPublishedbyWileyPeriodicals,Inc.
Key words: CAPRI;CASP; oligomerstate; blind prediction;protein interaction; proteindocking.
INTRODUCTION Computational approaches play a major role in all
|     |     |     |     |     |     |     | these endeavors. |     | Of particular |     | importance | are | methods |
| --- | --- | --- | --- | --- | --- | --- | ---------------- | --- | ------------- | --- | ---------- | --- | ------- |
Most cellular processes are carried out by physically for deriving accurate structural models of multiprotein
1
| interacting | proteins. | Characterizing |     | protein | interactions |     |             |          |      |            |             |     |        |
| ----------- | --------- | -------------- | --- | ------- | ------------ | --- | ----------- | -------- | ---- | ---------- | ----------- | --- | ------ |
|             |           |                |     |         |              |     | assemblies, | starting | from | the atomic | coordinates |     | of the |
and higher order assemblies is therefore a crucial step in individual components, the so-called “docking” algo-
| gaining an | understanding |     | of how | cells function. |     |     |         |         |            |           |          |     |          |
| ---------- | ------------- | --- | ------ | --------------- | --- | --- | ------- | ------- | ---------- | --------- | -------- | --- | -------- |
|            |               |     |        |                 |     |     | rithms, | and the | associated | energetic | criteria | for | singling |
Regrettably, protein assemblies are still poorly repre- 11–13
|     |     |     |     |     |     |     | out stable | binding | modes. |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | ---------- | ------- | ------ | --- | --- | --- | --- |
2
sented in the Protein Databank (PDB). Determining the Taking its inspiration from CASP, the community-
| structures | of such | assemblies | has | so far been | hampered | by  |                 |     |        |          |            |     |           |
| ---------- | ------- | ---------- | --- | ----------- | -------- | --- | --------------- | --- | ------ | -------- | ---------- | --- | --------- |
|            |         |            |     |             |          |     | wide initiative |     | on the | Critical | Assessment | of  | Predicted |
the difficulty in obtaining suitable crystals and diffraction Interactions (CAPRI), established in 2001, has been
data. But this limitation is being circumvented with the designed to test the performance of docking algorithms
| advent of | new powerful |     | electron | microscopy | techniques, |     |                                        |     |     |     |      |     |          |
| --------- | ------------ | --- | -------- | ---------- | ----------- | --- | -------------------------------------- | --- | --- | --- | ---- | --- | -------- |
|           |              |     |          |            |             |     | (http://www.ebi.ac.uk/msd-srv/capri/). |     |     |     | Just | as  | CASP has |
which now enable the structure determinations of very fostered the development of methods for the prediction
| large macromolecular |     | assemblies |     | at atomic | resolutions. | 3   |            |             |     |       |            |     |           |
| -------------------- | --- | ---------- | --- | --------- | ------------ | --- | ---------- | ----------- | --- | ----- | ---------- | --- | --------- |
|                      |     |            |     |           |              |     | of protein | structures, |     | CAPRI | has played | an  | important |
On the other hand, the repertoire of individual protein role in advancing the field of modeling protein assem-
3D structures has been increasingly filled, thanks to blies. Initially focusing on protein–protein docking and
| large-scale | structural | genomics |     | projects | such as | the PSI |         |             |       |     |          |     |            |
| ----------- | ---------- | -------- | --- | -------- | ------- | ------- | ------- | ----------- | ----- | --- | -------- | --- | ---------- |
|             |            |          |     |          |         |         | scoring | procedures, | CAPRI | has | expanded | its | horizon by |
(http://sbkb.org/) and others (http://www.thesgc.org/). including targets representing protein-peptide and pro-
| Given a | newly sequenced |     | protein, | the odds | are high | that |              |       |            |     |              |     |           |
| ------- | --------------- | --- | -------- | -------- | -------- | ---- | ------------ | ----- | ---------- | --- | ------------ | --- | --------- |
|         |                 |     |          |          |          |      | tein nucleic | acids | complexes. | It  | has moreover |     | conducted |
its 3D structure can be readily extrapolated from struc- experiments aimed at evaluating the ability of computa-
tures of related proteins deposited in the PDB. 4,5 More- tional methods to estimate binding affinity of protein–
14–16
over, thanks to the recent explosion of the number of protein complexes and to predict the positions of
17
available protein sequences, it is now becoming possible water molecules at the interfaces of protein complexes.
to model the structures of individual proteins with Considering the importance of macromolecular assem-
6,7
increasing accuracy from sequence information alone blies, and the new opportunities offered by the
as will be highlighted in the CASP11 results in this issue. recent progress in both experimental and computational
Structures from this increasingly rich repertoire may be techniques to probe and model these assemblies, a better
used as templates or scaffolds in protein design projects integration of the different computational approaches for
that have useful medical applications. 8,9 Larger protein modeling macromolecular assemblies and their building
assemblies can be modeled by integrating information on blocks was called for. Establishing closer ties between the
individual structures with various other types of data CASP and CAPRI communities appeared as an impor-
10
with the help of hybrid modeling techniques. tant step in this direction, inaugurated by running a
|     |     |     |     |     |     |     |     |     |     |     |     | PROTEINS | 325 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | -------- | --- |

M.F.Lensinketal.
joint CASP-CAPRI prediction experiment in the summer the targets (23) the goal was to model the interface (or
of 2014. The results of this experiment were summarized interfaces in the case of tetramers) between identical sub-
at the CASP11 meeting held in Dec 2014 in Cancun units, whose size varied between 44 and 669 residues but
Mexico, and are presented in detail in this report. was of (cid:2)250 residues on average. The majority of the
The CASP11-CAPRI experiment, representing CAPRI targets were obtained from structural genomics consortia.
Round 30, comprised 25 targets for which predictions of They represented mainly microbial proteins, whose func-
protein complexes were assessed. These targets repre- tion was often annotated as putative.
sented a subset of the 100 regular CASP11 targets. This Since it is not uncommon for docking approaches
subset comprised only “easy” CASP targets, those whose to use information on the symmetry of the complex to
3D structure could be readily modeled using standard restrain or filter docking poses, predictors needed to
homology modeling techniques. Targets that required be given reliable information on the biologically/func-
more sophisticated approaches (ab-initio modeling, or tionally relevant oligomeric state of the target complex
homology modeling using very distantly related tem- to be predicted. While self association between pro-
plates) were not considered, as the CAPRI community teins is common, with between 50 and 75% of pro-
had little experience with these approaches. The vast teins forming dimers in the cell, 20,21 this association
majority of the targets were homo-oligomers. CAPRI depends on the binding affinity between the subunits
groups were given the choice of modeling the subunit and on their concentration. Information on the oligo-
structures of these complexes themselves, or using mod- meric state is in principle derived using experimental
els made available by CASP participant, in time of the methods such as gel filtration or small-angle X-ray
docking calculations. scattering (SAXS), 22 and is usually communicated by
On average, about 25 CAPRI groups, and about 7 the authors upon submission of the atomic coordinates
CASP groups submitted docking predictions for each tar- to the PDB. With a majority of the targets being
get. About 12 CAPRI scorer groups per target partici- offered by structural genomics consortia before their
pated in the CAPRI scoring experiment, where coordinates were deposited in the PDB, author-
participants are invited to single out correct models from assigned oligomeric states were available to predictors
an ensemble of anonymized predicted complexes gener- only for a subset ((cid:2)15) of the targets, and those were
ated during the docking experiment. often tentative. For the remaining targets, the oligo-
In total, these groups submitted >9500 models that meric state was inferred from the crystal contacts using
were assessed against the 3D structures of the corre- the PISA software, 23 which although being a widely
sponding targets. The assessment was performed by the used standard in the field, may still yield erroneous
CAPRI assessment team, using the standard CAPRI assignments in a non-negligible fraction of the cases,
model quality measures.18,19 A major issue for the as will be shown in this analysis. Such incorrect
assessment, and for the Round as a whole, was the
assignments represented a confounding factor in this
uncertainties in the oligomeric state assignments for a
CAPRI round, but also allowed to show that docking
significant number of the targets. For many of these the
calculations may help to correct them.
assigned state at the time of the experiment was inferred
solely from the crystal contacts by computational meth-
ods, which can be unreliable.
GLOBAL OVERVIEW OF THE
In presenting the CAPRI Round 30 assessment results
PREDICTION EXPERIMENT
here, we highlight this issue and the more general chal-
lenge of correctly predicting the association modes of As in typical CAPRI Rounds, CAPRI predictor groups
weaker complexes of identical subunits, and those of were provided with the amino-acid sequence of the tar-
higher order homo-oligomers. In addition, we examine get protein (for homo-oligomers), or proteins (for heter-
the influence of the accuracy of the modeled subunits on ocomplexes), and with some relevant details about the
the performance of the docking and scoring predictions, protein, communicated by the structural biologists. Using
and evaluate the extent to which docking procedures the sequence information, the groups were then invited
confer an advantage over standard homology modeling to model the 3D structure of the protein or proteins,
methods in predicting homo-oligomer complexes. and to derive the atomic structure of the complex. To
help with the homology-modeling task, with which
CASP participants are usually more experienced than
THE TARGETS
their CAPRI colleagues, 3D models of individual target
The 25 targets of the joint CASP-CAPRI experiment proteins predicted by CASP participants were made
T1 are listed in Table I. Of these 23 are homo-oligomers, available to CAPRI groups for use in their docking calcu-
with 18 declared to be dimers and five to be tetramers, lations. A good number of CAPRI groups, but not all,
and two heterocomplexes. Hence for the majority of took up this offer.
326 PROTEINS

PredictionofHomoandHeteroproteinComplexesbyProteinDockingandModeling
TableI
TheCAPRI-CASP11TargetsofCAPRIRound30
| TargetID |     |     | Quaternarystate |     |     |     |     |     |
| -------- | --- | --- | --------------- | --- | --- | --- | --- | --- |
Buried
area((cid:2)2)
| CAPRI | CASP | Contributor | Author | PISA | Residues |     |     | Protein |
| ----- | ---- | ----------- | ------ | ---- | -------- | --- | --- | ------- |
T68 T0759 NSGC 1or2 1 109 860 Plectin1and2repeats(HR9083A)ofthe
HumanPeriplakin
T69 T0764 JCSG 2 2 341 2415 Putativeesterase(BDI_1566)fromPara-
bacteroidesdistasonis
T70 T0765 JCSG 2 4 128 2030 ModulatorproteinMzrA(KPN_03524)from
Klebsiellapneumoniaesubsp.
| T71 | T0768 | JCSG | 4   | 4   | 170 | 2380 | Leucinerichrepeatprotein(BAC- |     |
| --- | ----- | ---- | --- | --- | --- | ---- | ----------------------------- | --- |
CAP_00569)fromBacteroidescapillosus
ATCC29799
T72 T0770 JCSG 2 2 488 1120 SusDhomolog(BT2259)fromBacteroides
thetaiotaomicron
T73 T0772 JCSG 4 4 265 5900 Putativeglycosylhydrolase(BDI_3914)
fromParabacteroidesdistasonis
T74 T0774 JCSG 1 4 379 2040 Hypotheticalprotein(BVU_2522)from
Bacteroidesvulgatus
T75 T0776 JCSG 2 2 256 1040 PutativeGDSL-likelipase(PARMER_00689)
fromParabacteroidesmerdae(ATCC
43184)
T77 T0780 JCSG 2 2 259 1600 Conservedhypotheticalprotein(SP_1560)
fromStreptococcuspneumoniaeTIGR4
T78 T0786 Non-SGI 4 4 264 4160 Hypotheticalprotein(BCE0241)fromBacil-
luscereus
| T79 | T0792 | Non-SGI |     | 2   | 80  | 680 | OSKAR-N |     |
| --- | ----- | ------- | --- | --- | --- | --- | ------- | --- |
T80 T0801 NPPB 2 2 376 1960 SugaraminotransferaseWecEfromEsch-
erichiacoliK-12
T81 T0797 Non-SGI 2 2 44 1070 cGMP-dependentproteinkinaseIIleucine
zipper
|     | T0798 |         | 2   | 2   | 198 |      | Rab11bprotein            |     |
| --- | ----- | ------- | --- | --- | --- | ---- | ------------------------ | --- |
| T82 | T0805 | Non-SGI | 2   | 2   | 214 | 3250 | Nitro-reductaserv3368    |     |
| T84 | T0811 | NYSGRC  |     | 2   | 255 | 1740 | Triosephosphateisomerase |     |
T85 T0813 NYSGRC 2 2 307 4620 Cyclohexadienyldehydrogenasefrom
Sinorhizobiummelilotiincomplexwith
NADP
T86 T0815 NYSGRC 2 2 106 470 Putativepolyketidecyclase(protein
SMa1630)fromSinorhizobiummeliloti
T87 T0819 NYSGRC 2 2 373 3430 Histidinol-phosphateaminotransferase
fromSinorhizobiummelilotiincomplex
withpyridoxal-5'-phosphate
| T88 | T0825 | Non-SGI | 2   | 2   | 205 | 1350 | WRAP-5                              |     |
| --- | ----- | ------- | --- | --- | --- | ---- | ----------------------------------- | --- |
| T89 | T0840 | Non-SGI | 1   |     | 669 | 870  | RONreceptortyrosinekinasesubunit    |     |
|     | T0841 |         | 1   |     | 253 |      | Macrophagestimulatingproteinsubunit |     |
(MSP)
| T90 | T0843 | MCSG | 2   | 2   | 369 | 2360 | Ats13         |     |
| --- | ----- | ---- | --- | --- | --- | ---- | ------------- | --- |
| T91 | T0847 | SGC  | 1   | 2   | 176 | 1320 | HumanBj-Tsa-9 |     |
T92 T0849 MCSG 2 2 240 1900 GlutathioneS-transferasedomainfrom
HaliangiumochraceumDSM14365
T93 T0851 MCSG 2 2 456 2680 Cals8fromMicromonosporaechinospora
(P294Smutant)
| T94 | T0852 | MCSG | 2   | 2   | 414 | 1190 | APC103154 |     |
| --- | ----- | ---- | --- | --- | --- | ---- | --------- | --- |
BoldnumbersunderQuaternaryStateindicatetheoligomericstateassignmentsavailableatthetimeofthepredictionexperiment;1(monomer),2(dimer),4(tet-
ramer);numbersinregularfontsindicatesubsequentassignmentscollectedfromthePDBentriesforthetargetstructures.
NSGC, Northeast Structural Genomics Consortium; JCSG, Joint Center for Structural Genomics; Non-SGI, Non-SGI research Centers and others; NNPB, NatPro
PSI:Biology;NYSGRC,NewYorkStructuralGenomicsResearchCenter;MCSG,MidwestCenterforStructuralGenomics;SGC,StructuralGenomicsConsortium.
In addition to submitting 10 models for each target the ensemble of uploaded models using the scoring
complex, predictors were invited to upload a set of function of their choice, and submit their own 10
100 models. Once all the submissions were completed, best ranking ones. The typical timelines per target
the uploaded models were shuffled and made available were about 3 weeks for the homology modeling and
to all groups as part of the CAPRI scoring experiment. docking predictions, and 3 days for the scoring
| The “scorer” | groups | were in | turn invited | to evaluate | experiment. |     |     |     |
| ------------ | ------ | ------- | ------------ | ----------- | ----------- | --- | --- | --- |
PROTEINS 327

M.F.Lensinketal.
TableII
CAPRIRound30ExperimentStatistics
|     |          |     |     |     |     | Numberofgroups |     |     |      |     | Numberofmodels |     |      |
| --- | -------- | --- | --- | --- | --- | -------------- | --- | --- | ---- | --- | -------------- | --- | ---- |
|     | TargetID |     |     |     |     | CAPRI          |     |     |      |     | CAPRI          |     |      |
|     |          |     |     |     |     |                |     |     | CASP |     |                |     | CASP |
a
CAPRI CASP PDB Predictors Uploaders Scorers Predictors Predictors Uploaders Scorers Predictors
| T68 |     | T0759 | 4q28 | 2                     | 23  | 10  |     | 12  | 3   | 221 | 1000 | 120 | 7   |
| --- | --- | ----- | ---- | --------------------- | --- | --- | --- | --- | --- | --- | ---- | --- | --- |
| T69 |     | T0764 | 4q34 | 2                     | 28  | 10  |     | 14  | 7   | 266 | 1000 | 132 | 17  |
| T70 |     | T0765 | 4pwu | 2                     | 23  | 8   |     | 13  | 5   | 221 | 710  | 130 | 18  |
| T71 |     | T0768 | 4oju | 3                     | 22  | 9   |     | 14  | 1   | 214 | 810  | 131 | 1   |
| T72 |     | T0770 | 4q69 | 3                     | 25  | 11  |     | 13  | 4   | 244 | 914  | 130 | 11  |
| T73 |     | T0772 | 4qhz | 2                     | 23  | 11  |     | 11  | 7   | 221 | 1195 | 110 | 16  |
| T74 |     | T0774 | 4qb7 | 2                     | 22  | 11  |     | 10  | 7   | 202 | 911  | 96  | 11  |
| T75 |     | T0776 | 4q9a | 1                     | 26  | 12  |     | 12  | 8   | 253 | 840  | 120 | 21  |
| T76 |     | T0779 |      | Cancelled–nostructure |     |     |     |     |     |     |      |     |     |
| T77 |     | T0780 | 4qdy | 4                     | 24  | 12  |     | 12  | 6   | 229 | 971  | 120 | 12  |
| T78 |     | T0786 | 4qvu | 2                     | 24  | 10  |     | 11  | 5   | 229 | 818  | 110 | 15  |
| T79 |     | T0792 | 5a49 | 3                     | 25  | 11  |     | 12  | 9   | 242 | 900  | 120 | 23  |
| T80 |     | T0801 | 4piw | 1                     | 27  | 10  |     | 12  | 8   | 264 | 911  | 120 | 27  |
| T81 |     | T0797 | 4ojk | 1                     | 23  | 9   |     | 11  | 20  | 218 | 641  | 110 | 64  |
T0798
| T82 |     | T0805 | b   | 1                                         | 25  | 10  |     | 12  | 9   | 242 | 911 | 120 | 27  |
| --- | --- | ----- | --- | ----------------------------------------- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| T83 |     | T0809 |     | Cancelled–articlefromdifferentgrouponline |     |     |     |     |     |     |     |     |     |
b
| T84 |     | T0811 |      | 1   | 25  | 10  |     | 12  | 10  | 241 | 910  | 120 | 28  |
| --- | --- | ----- | ---- | --- | --- | --- | --- | --- | --- | --- | ---- | --- | --- |
| T85 |     | T0813 | 4wji | 1   | 25  | 11  |     | 12  | 8   | 241 | 920  | 120 | 21  |
| T86 |     | T0815 | 4u13 | 2   | 26  | 11  |     | 12  | 9   | 251 | 1010 | 119 | 25  |
| T87 |     | T0819 | 4wbt | 1   | 24  | 10  |     | 12  | 9   | 231 | 894  | 120 | 25  |
| T88 |     | T0825 | b    | 1   | 27  | 10  |     | 13  | 18  | 261 | 910  | 130 | 62  |
| T89 |     | T0840 | b    | 1   | 22  | 9   |     | 11  | 55  | 211 | 790  | 110 | 243 |
T0841
| T90 |     | T0843 | 4xau | 1   | 23  | 9   |     | 11  | 9   | 221 | 811 | 110 | 28  |
| --- | --- | ----- | ---- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| T91 |     | T0847 | 4urj | 1   | 25  | 9   |     | 11  | 9   | 242 | 798 | 110 | 24  |
| T92 |     | T0849 | 4w66 | 1   | 23  | 9   |     | 11  | 9   | 225 | 789 | 110 | 33  |
| T93 |     | T0851 | 4wb1 | 1   | 22  | 9   |     | 11  | 8   | 213 | 697 | 110 | 27  |
| T94 |     | T0852 | 4w9r | 1   | 22  | 9   |     | 12  | 8   | 215 | 783 | 120 | 21  |
Thenumberofgroupscorrespondstoregisteredgroupsthateffectivelysubmittedmodelsfortherespectivetarget.Thenumberofmodelsrepresentssubmittedmodels,
regardlessofqualityandincludesdisqualifiedmodels.CAPRIgroupsareallowedtosubmitnomorethantheir10bestmodels,whereasCASPgroupsareallowedto
submitnomorethantheir5bestmodels.
aNumberofinterfacesassessed.
bNotyetreleased.
T2 Table II lists for each target the number of groups sub- to groups solely interested in testing their scoring
| mitting    | predictions |           | and     | the number |            | of models | assessed.  | functions. |     |     |     |     |     |
| ---------- | ----------- | --------- | ------- | ---------- | ---------- | --------- | ---------- | ---------- | --- | --- | --- | --- | --- |
| On         | average     | (cid:2)25 | CAPRI   | groups     | submitted  |           | a total of |            |     |     |     |     |     |
| (cid:2)230 | models      | per       | target, | and        | an average | of        | 12 scorer  |            |     |     |     |     |     |
groups submitted a total of (cid:2)120 models per target. SYNOPSIS OF THE PREDICTION
| With   | the        | exception | of    | three targets, | an        | average    | of seven | METHODS |                       |     |             |          |          |
| ------ | ---------- | --------- | ----- | -------------- | --------- | ---------- | -------- | ------- | --------------------- | --- | ----------- | -------- | -------- |
| groups | registered |           | with  | CASP           | submitted | a total    | of any-  |         |                       |     |             |          |          |
|        |            |           |       |                |           |            |          |         | Round 30 participants |     | used a wide | range of | modeling |
| where  | between    |           | 1 and | 33 models      | for       | individual | targets. |         |                       |     |             |          |          |
CASP predictors participated in larger numbers in the methods and software tools to generate the submitted
|            |     |        |         |     |        |                 |     | models. | In addition, |     | the approaches | used by | a given |
| ---------- | --- | ------ | ------- | --- | ------ | --------------- | --- | ------- | ------------ | --- | -------------- | ------- | ------- |
| prediction |     | of T88 | (T0825) | and | of the | heterocomplexes |     |         |              |     |                |         |         |
(T89 – T0840/T0841 and T81 – T0797/T0798), group often differed across targets. Here, we provide
where the CASP targets were defined as the oligomeric only a short synopsis of the main methodological
structures. approaches. For a more detailed description of the meth-
Table II also lists the uploader groups and the models ods and modeling strategies, readers are referred to the
|      |      |      |           |         |         |            |      | extended | Methods | Abstracts | provided | by individual | par- |
| ---- | ---- | ---- | --------- | ------- | ------- | ---------- | ---- | -------- | ------- | --------- | -------- | ------------- | ---- |
| that | they | make | available | for the | scoring | experiment | (100 |          |         |           |          |               |      |
models per target per uploader group). As detailed ticipants (see Supporting Information Table S6).
above, the uploaded models are complexes output by the Templates, representing known structures of homologs
docking calculations carried out by individual partici- to a given target, stored in the PDB, were used in a
pants for a given target. Models, uploaded by the differ- number of ways. Most commonly, they were employed
ent groups, are anonymized, shuffled, and made available to model the 3D structures of individual subunits. Some
328 PROTEINS

PredictionofHomoandHeteroproteinComplexesbyProteinDockingandModeling
C
O
L
O
R
Figure 1
SchematicillustrationoftheCAPRIassessmentcriteria.Thefollowingquantitieswerecomputedforeachtarget:(1)alltheresidue-residuecontacts
betweentheReceptor(R)andtheLigand(L),and(2)theresiduescontributingtotheinterfaceofeachofthecomponentsofthecomplex.Inter-
faceresiduesweredefinedonthebasisoftheircontributiontotheinterfacearea,asdescribedinreferences. 18,19 Foreachsubmittedmodelthefol-
lowingquantitieswerecomputed:thefractionsf(nat)ofnativeandf(non-nat)ofnon-nativecontactsinthepredictedinterface;therootmean
squaredisplacement(rmsd)ofthebackboneatomsoftheligand(L-rms),themis-orientationangleh andtheresidualdisplacementd ofthe
|     |     |     |     |     |     |     |     | L   |     | L   |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
ligandcenterofmass,afterthereceptorinthemodelandexperimentalstructureswereoptimallysuperimposed.InadditionwecomputedI-rms,
thermsdofthebackboneatomsofallinterfaceresiduesaftertheyhavebeenoptimallysuperimposed.Heretheinterfaceresiduesweredefinedless
stringentlyonthebasisofresidue-residuecontacts(seeRefs.18,19).
CAPRI participants selected their own templates and guide the docking calculations or to select docking solu-
used a variety of custom built or well-established algo- tions. Others used the dimeric templates directly to model
|     |     |     |     | 24  |     | 25  |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
rithms such as Modeller, Swiss-Model, or the target dimer (template-based “docking” 32–34 ). Less
| ROSETTA, | 26 to | model | the | subunit | structures. | Others |                  |                    |                |      |
| -------- | ----- | ----- | --- | ------- | ----------- | ------ | ---------------- | ------------------ | -------------- | ---- |
|          |       |       |     |         |             |        | than a hand-full | of groups employed | template-based | mod- |
used the models produced by various servers participat- eling alone for all or most of the targets.
ing in the CASP11 experiment and made available to To model tetrameric targets, most groups proceeded in
CAPRI groups, or servers of other groups (HAD- two steps. They used either known dimeric homologs, or
27
DOCK ). The quality of the CASP server models was docking methods to build the dimer portion of the tet-
| usually first | assessed |     | using various | criteria | and | the best |                 |                   |            |             |
| ------------- | -------- | --- | ------------- | -------- | --- | -------- | --------------- | ----------------- | ---------- | ----------- |
|               |          |     |               |          |     |          | ramer, and then | run their docking | procedures | to generate |
quality models were selected for the docking calculations. a dimer-of-dimers, representing the predicted tetramer.
| Some groups   | selected  |         | a single     | best model   | for a       | given tar- |              |            |     |     |
| ------------- | --------- | ------- | ------------ | ------------ | ----------- | ---------- | ------------ | ---------- | --- | --- |
| get, whereas  | others    | used    | several      | models       | (sometimes  | up         |              |            |     |     |
|               |           |         |              |              |             |            | ASSESSMENT   | PROCEDURES |     |     |
| to five       | models).  | Several | groups       | additionally |             | used loop  |              |            |     |     |
|               |           |         |              |              |             |            | AND CRITERIA |            |     |     |
| modeling      | to adjust | the     | conformation |              | of loops    | regions,   |              |            |     |     |
| and subjected | the       | subunit | models       | to energy    | refinement. |            |              |            |     |     |
ThestandardCAPRIassessmentprotocol
| The majority |             | of      | CAPRI | participants | used | protein   |               |          |                 |      |
| ------------ | ----------- | ------- | ----- | ------------ | ---- | --------- | ------------- | -------- | --------------- | ---- |
|              |             |         |       |              |      |           | The predicted | homo and | heterocomplexes | were |
| docking      | and scoring | methods |       | to generate  | and  | rank can- |               |          |                 |      |
didate complexes. Many employed their own docking assessed by the CAPRI assessment team, using the stand-
methods, some of which were designed to handle sym- ard CAPRI assessment protocol, which evaluates the cor-
metric assemblies, whereas others relied on well- respondence between predicted complex and the target
18,19
| established | docking | algorithms |     | such as | HEX, 28 | ZDock, 29 | structure. |     |     |     |
| ----------- | ------- | ---------- | --- | ------- | ------- | --------- | ---------- | --- | --- | --- |
30
RosettaDock, as well as on docking programs such as This protocol (summarized in Fig. 1) first defines the F1
MZDock31
which apply symmetry constraints. set of residues common to all the submitted models and
When templates were available for a given target the target, so as to enable the comparison of residue-
(mostly for homodimers), some participants used the dependent quantities, such as the root mean square devi-
information from these templates (consensus interface res- ation (rmsd) of the models versus the target structure.
idues, contacts, or relative arrangement of subunits) to Models where the sequence identity to the target is too
|     |     |     |     |     |     |     |     |     | PROTEINS | 329 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | -------- | --- |

M.F.Lensinketal.
TableIII inconsistent (Table I). Only about 15 targets had an oli-
SummaryofCAPRICriteriaforRankingPredictedComplexes gomeric state assigned by the authors at the time of the
experiment.
Score f(nat) L-rms I-rms
To address this problem in the assessment, the PISA
*** High (cid:3)0.5 (cid:4)1.0 OR (cid:4)1.0
software was used to generate all the crystal contacts for
** Medium (cid:3)0.3 <1.0–5.0] OR <1.0–2.0]
* Acceptable (cid:3)0.1 <5.0–10.0] OR <2.0–4.0] each target and to compute the corresponding interface
Incorrect <0.1 >10.0 AND >4.0 areas. The interfaces were then ranked according to size
of the interface. In candidate dimer targets, submitted
models were usually evaluated against 1 or 2 of the larg-
low are not assessed. The threshold is determined on a
est interfaces of the target, and acceptable or better mod-
per-target basis, but is typically set to 70%.
els for any or all of these interfaces were tallied. For
The set of common residues is used to evaluate the
candidate tetramer targets, the relevant largest interfaces
two main rmsd-based quantities used in the assessment:
for each target were identified in the crystal structure,
the ligand rmsd (L-rms) and the interface rmsd (I-rms).
and predicted models were evaluated by comparing in
L-rms is the backbone rmsd over the common set of
turn each pair of interacting subunits in the model to
ligand residues after a structural superposition of the
each of the relevant pairs of interacting subunits in the
receptor. I-rms is the backbone rmsd calculated over the
target (Supporting Information Fig. S1), and again the
common set of interface residues after a structural super-
best predicted interfaces were retained for the tally. One
position of these residues. An interface residue is defined
of the two bonafide heterocomplexes was also evaluated
as such when any of its atoms (hydrogens excluded) are
against multiple interfaces.
found within 10 A˚ of any of the atoms of the binding
partner.
Evaluatingtheaccuracyofthe3Dmodelsof
An important third quantity whereby models are individualsubunits
assessed is f(nat), representing the fraction of native con-
Since this experiment was a close collaboration
tacts in the target, that is, reproduced in the model. This
between CAPRI and CASP, the quality of the 3D models
quantity takes all the protein residues into account. A
of individual subunits in the predicted complexes was
ligand-receptor contact is defined as any pair of ligand-
35
receptor atoms within 5 A˚ distance. Atomic contacts assessed by the CASP team using the LGA program,
below 3 A˚ are considered as clashes; predictions with too which is the basic tool for model/target comparison in
36,37
CASP. The tool can be run in two evaluation
many clashes are disqualified. The clash threshold varies
modes. In the sequence-dependent mode, the algorithm
with the target and is defined as the average number of
assumes that each residue in the model corresponds to a
clashes in the set of predictions plus two standard devia-
residue with the same number in the target, while in the
tions. The quantities f(nat), L-rms and I-rms together
sequence-independent mode this restriction is not
determine the quality of a predicted model, and based
applied. The program searches for optimal superimposi-
on those three parameters models are ranked into four
tions between two structures at different distance cutoffs
categories: High quality, medium quality, acceptable
and returns two main accuracy scores; GDT_TS and
T3 quality and incorrect, as summarized in Table III.
LGA_S. The GDT_TS score is calculated in the sequence-
dependent mode and represents the average percentage
ApplyingtheCAPRIassessmentprotocolto
of residues that are in close proximity in two structures
homo-oligomers
optimally superimposed using four selected distance cut-
Evaluating models of homo and heteroprotein com- offs (see Ref. 38 for details). The LGA_S score is calcu-
plexes against the corresponding target structure is a lated in both evaluation modes and represents a
well-defined problem when the target complex is unam- weighted sum of the auxiliary LCS and GDTscores from
biguously defined, for example, if the target association the superimpositions built for the full set of distance cut-
mode and corresponding interface represents the biologi- offs (see Ref. 35 for details). We have run the evaluation
cally relevant unit. This is usually, although not always, in both modes, but since the CAPRI submission format
the case for binary heterocomplexes, but was not the sit- permits different residue numbering, we used the LGA_S
uation encountered in this experiment for the homo- score from the sequence-independent analysis as the
oligomer targets. All except two of the 25 targets for main measure of the subunit accuracy assessment. This
which predictions were evaluated here represent homo- score is expressed on a scale from 0 to 100, with 100 rep-
oligomers. For about half of these targets the oligomeric resenting a model that perfectly fits the target. The rmsd
state was deemed unreliable, as it was either only values for subunit models cited throughout the text are
inferred computationally from the crystal structure using those computed by LGA software. We verified that for
23
the PISA software or because the authors’ assignment about 80% of the assessed models the GDT-TS and
and inferred oligomeric states, although available, were LGA-S scores differed by <15 units, indicating that these
330 PROTEINS

PredictionofHomoandHeteroproteinComplexesbyProteinDockingandModeling
models correspond to near identical structural align- the experimental oligomer structure: A similar procedure
ments with the corresponding targets, in line with the to those described above was applied. Although this
fact that the majority of the targets of this Round repre- time, the best templates were identified by searching for
sent proteins that could be readily modeled by homol- proteins with the highest structural similarity to the tar-
ogy. Of the remaining 20% with larger differences get oligomer structure. The search was performed using
44
between the 2 scores, 18% correspond to disqualified the multimeric structure alignment tool MM-align.
models or incorrect complexes and 2% correspond to For computational efficiency, MM-align was applied only
acceptable (or higher quality) predicted complexes. Their to the 100 proteins with the highest monomer structure
impact on the analysis is therefore negligible. similarity to the target. Models were built using MOD-
|     |     |     |     |     |     |     | ELLER based | on  | the alignment |     | output | by MM-align. |     |     |
| --- | --- | --- | --- | --- | --- | --- | ----------- | --- | ------------- | --- | ------ | ------------ | --- | --- |
Buildingtargetmodelsbasedonthebest
availabletemplates
RESULTS
| In order | to better | estimate |     | the added | value | of protein |     |     |     |     |     |     |     |     |
| -------- | --------- | -------- | --- | --------- | ----- | ---------- | --- | --- | --- | --- | --- | --- | --- | --- |
docking procedures and template-based modeling techni- This section is divided into three parts. The first part
|         |        |             |     |       |            |         | presents | the prediction |     | results | for the | 25  | individual | tar- |
| ------- | ------ | ----------- | --- | ----- | ---------- | ------- | -------- | -------------- | --- | ------- | ------- | --- | ---------- | ---- |
| ques it | seemed | of interest | to  | build | a baseline | against |          |                |     |         |         |     |            |      |
which the different approaches could be benchmarked. gets for which the docking and scoring experiments were
To this end, the best oligomeric structure template for conducted. In the second part, we present an overview of
each target available at the time of the predictions was the results across targets and across predictor and scorer
identified. Based on this template, the target model was groups, respectively. In the third part, we review the
|             |            |          |     |            |     |           | accuracy | of the models |     | of individual |     | subunits | in  | the pre- |
| ----------- | ---------- | -------- | --- | ---------- | --- | --------- | -------- | ------------- | --- | ------------- | --- | -------- | --- | -------- |
| built using | a standard | modeling |     | procedure, | and | the qual- |          |               |     |               |     |          |     |          |
ity of this model was assessed using the CAPRI evalua- dicted oligomers, and how this accuracy influences the
|               |           |            |     |         |           |       | performance | of docking |     | procedures. |     |     |     |     |
| ------------- | --------- | ---------- | --- | ------- | --------- | ----- | ----------- | ---------- | --- | ----------- | --- | --- | --- | --- |
| tion criteria | described | above.     |     |         |           |       |             |            |     |             |     |     |     |     |
| To identify   | the       | templates, | the | protein | structure | data- |             |            |     |             |     |     |     |     |
base “PDB70” containing proteins of mutual sequence Predictionresultsforindividualtargets
39
| identity | (cid:4)70%  | was downloaded |        | from | HHsuite.   | The  |                |     |          |      |      |           |      |      |
| -------- | ----------- | -------------- | ------ | ---- | ---------- | ---- | -------------- | --- | -------- | ---- | ---- | --------- | ---- | ---- |
|          |             |                |        |      |            |      | Easy homodimer |     | targets: | T69, | T75, | T80, T82, | T84, | T85, |
| database | was updated | twice          | during | the  | experiment | (See |                |     |          |      |      |           |      |      |
T87,T90,T91,T92,T93,T94
| Supporting | Information |     | Table | S5 for | the release | date of |     |     |     |     |     |     |     |     |
| ---------- | ----------- | --- | ----- | ------ | ----------- | ------- | --- | --- | --- | --- | --- | --- | --- | --- |
the database used for each target). Only homo-complexes The 12 targets in this category comprised some of the
|                 |           |           |           |               |     |               | largest subunits |         | of the | entire  | evaluated | target    | set, | with   |
| --------------- | --------- | --------- | --------- | ------------- | --- | ------------- | ---------------- | ------- | ------ | ------- | --------- | --------- | ---- | ------ |
| were considered |           | for this  | analysis. |               |     |               |                  |         |        |         |           |           |      |        |
|                 |           |           |           |               |     |               | sizes ranging    | between |        | 176 and | 456       | residues. | Four | of the |
| The best        | available | templates |           | were detected |     | in three dif- |                  |         |        |         |           |           |      |        |
ferent ways and target models were generated from the targets were multi-domain proteins (T85, T87, T90, and
|           |             |     |               |     |       |             | T93), and | one (T82) | was | an intertwined |     | dimer. |     |     |
| --------- | ----------- | --- | ------------- | --- | ----- | ----------- | --------- | --------- | --- | -------------- | --- | ------ | --- | --- |
| templates | as follows: |     | (1) Detection |     | based | on sequence |           |           |     |                |     |        |     |     |
information alone: For each target sequence, proteins In the following, we present examples of the perform-
|            |            |      |          |     |        |             | ance achieved | for | this | category | of targets. | Detailed |     | results |
| ---------- | ---------- | ---- | -------- | --- | ------ | ----------- | ------------- | --- | ---- | -------- | ----------- | -------- | --- | ------- |
| related to | the target | were | searched |     | for in | the protein |               |     |      |          |             |          |     |         |
structure database by HHsearch40 in the local alignment for all the targets of Round 30 can be found in the Sup-
|           |             |     |            | 41    |     |             | porting | Information | Table | S2, | and on | the | CAPRI | website |
| --------- | ----------- | --- | ---------- | ----- | --- | ----------- | ------- | ----------- | ----- | --- | ------ | --- | ----- | ------- |
| mode with | the Viterbi |     | algorithm. | Among |     | the top 100 |         |             |       |     |        |     |       |         |
entries, up to 10 proteins that are in the desired (URL: http://www.ebi.ac.uk/msd-srv/capri/).
oligomer state were selected as templates. When more An illustrative example of the average performance
|          |          |            |     |      |           |            | obtained | for this | category | of  | targets | is that | obtained | for |
| -------- | -------- | ---------- | --- | ---- | --------- | ---------- | -------- | -------- | -------- | --- | ------- | ------- | -------- | --- |
| than two | assembly | structures |     | with | different | interfaces |          |          |          |     |         |         |          |     |
were identified, the best ranking one was selected as tem- target T69 (T0764): a 341-residue putative esterase
|            |        |     |          |           |      |         | (BDI_1566) | from | Parabacteroides |     | distasonis. |     | The | submit- |
| ---------- | ------ | --- | -------- | --------- | ---- | ------- | ---------- | ---- | --------------- | --- | ----------- | --- | --- | ------- |
| plate. The | target | and | template | sequences | were | aligned |            |      |                 |     |             |     |     |         |
using HHalign 40 in the global alignment mode with the ted models for this target were evaluated against two
|         |          |            |     |       |        |          | interfaces | in the | crystal | structure | of  | this protein, |     | gener- |
| ------- | -------- | ---------- | --- | ----- | ------ | -------- | ---------- | ------ | ------- | --------- | --- | ------------- | --- | ------ |
| maximum | accuracy | algorithm. |     | Based | on the | sequence |            |        |         |           |     |               |     |        |
alignments, oligomer models were built using MODEL- ated by applying the crystallographic symmetry opera-
LER. 42 The model with the lowest MODELLER energy tions listed in the Supporting Information Table S1, and
A˚2)
out of 10 models was selected for further analysis. (2) depicted in Figure 2(a): one large interface (2415 F2
Detection based on the experimental monomer structure: and a smaller interface (622 A˚2). Good prediction results
Proteins with highest structural similarity to the experi- were obtained only for interface 1. Twenty-eight CAPRI
mental monomer structure were searched for using TM- predictor groups submitted a total of 266 models for this
| 43     |       |         |     |          |       |             | homodimer. | Of  | these, | 30 were | of acceptable |     | quality | and |
| ------ | ----- | ------- | --- | -------- | ----- | ----------- | ---------- | --- | ------ | ------- | ------------- | --- | ------- | --- |
| align. | Among | the top | 100 | entries, | up to | 10 proteins |            |     |        |         |               |     |         |     |
that are in the desired oligomer state were selected as 57 were of medium quality. Twelve predictor groups and
templates as described above. Based on the target- three docking servers submitted at least one model of
template alignments output by TM-align, models were acceptable quality or better. Among those, nine groups
built using MODELLER, and the lowest energy model and one server (CLUSPRO) submitted at least 1 medium
was selected as described above. (3) Detection based on quality model. The best performance (10 medium quality
|     |     |     |     |     |     |     |     |     |     |     |     | PROTEINS |     | 331 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | -------- | --- | --- |

M.F.Lensinketal.
C
O
L
O
R
Figure 2
Targetstructuresandpredictionresultsforeasydimertargets.T69(T0764),aPutativeesterase(BDI_1566)fromParabacteroidesdistasonis,PDBcode
4Q34.(a)Targetstructure,withhighlightedinterfaces(1,2).(b)Globaldockingpredictionresultsdisplayingonesubunitincartoonrepresentation,
withthecenterofmassofthesecondsubunitinthetarget(redsphere),andindockingsolutionssubmittedbyCAPRIpredictors(lightbluespheres),
CAPRIscorers(darkbluespheres),andCASPpredictors(yellowspheres).T80(T0801),asugaraminotransferaseWecEfromEscherichiacoliK-12,
PDBcode4PIW.(c)Targetstructure.(d)Globaldockingpredictionresultsbydifferentpredictorgroups(seelegend(b)fordetail).T82(T0805)
Nitroreductase(structuresunreleased).(e)Targetstructure.(f)Globaldockingpredictionresultsbydifferentpredictorgroups.T94(T0852),unchar-
acterizedproteinCoch_1243fromCapnocytophagaochraceaDSM7271,PDBcode4W9R.(g)Targetstructure.(h)Globaldockingpredictionresults
bydifferentpredictorgroups.
models) was obtained by the groups of Seok, Lee and experiment, and thus not singling out even their own
Guerois, followed closely by the groups of Zou, Shen, best models from the uploaded anonymized set of pre-
and Eisenstein (see Supporting Information Table S2 for dicted complexes, highlighting yet again the distinct
the complete ranking) nature of the docking and scoring procedures.
The best model for this target, obtained by Guerois, An important factor in the successful predictions was
had an f(nat) value of 49%, and L-rms and I-rms values the overall good accuracy of the 3D models used by pre-
of 2.88 and 2.12 A˚, respectively (Supporting Information dictors in the docking calculations (see Fig. 6 and CAPRI
Table S4). website for detailed values). The best models had an
Six groups, registered with CASP, submitted in total LGA_S score of (cid:2)85 (backbone rmsd of (cid:2)3.9 A˚), and
12 models for this target, comprising one acceptable only a few models had LGA_S scores lower than 40
model by the group of Umeyama and one medium qual- (backbone rmsd>10 A˚) (values for all models are avail-
ity model by the Baker group. The global landscape of able on the CAPRI website). The accuracy of the 3D
all the predicted models by the different groups is out- models across targets and its influence on the predictions
lined in Figure 2(b). will be discussed in a dedicated section below.
An even better performance was achieved by the Very good predictions were obtained for T82 (T0805),
CAPRI scoring experiment (Supporting Information the nitro-reductase rv3368, a significantly intertwined
Table S2). Of the 14 groups participating in this experi- dimer with unstructured arms reaching out to the neigh-
ment, 12 submitted at least two models of medium qual- boring subunit and a subunit interface area of 3250 A˚2
ity. The best performance was achieved by Kihara (10 [Fig. 2(e,f)]. The majority of the models of the individ-
medium quality models), closely followed by Zou and ual subunits were quite accurate with LGA_S values of
Grudinin, with eight and five medium quality models, 60–85 (backbone rmsd <5 A˚) (see CAPRI website). As
respectively. As already observed in previous CAPRI eval- many as 54 medium quality models and 17 acceptable
uations the best performers in the docking calculations models were submitted by CAPRI participants, 99 mod-
were not necessarily performing as well in the scoring els of acceptable quality or better were submitted by
332 PROTEINS

PredictionofHomoandHeteroproteinComplexesbyProteinDockingandModeling
CAPRI scorer groups, and 11 acceptable models or better monomer by the authors. The good docking perform-
were submitted by three CASP groups (Supporting Infor- ance for this target and the fact that the dimer interface
mation Table S2). The high success rate for both com- (1320 A˚2) is within the range expected for proteins of
plex predictions and subunit modeling stems from the this size (176 residues), 45 suggests that this protein
fact that most predictors made good use of known struc- forms a dimer.
tures of related homodimers in the PDB in which the A somewhat lower performance was achieved for T92
intertwining mode was well conserved. These known (T0849) the glutathione S-transferase domain from Hal-
dimer structures were mainly used in templates for mod- iangium ochraceum), and for T94 (T0852), an uncharac-
eling the target dimer (template-based docking). terized 2-domain protein (putative esterase according to
Very similar participation, number of submitted mod- Pfam) Coch_1243 from Capnocytophaga ochracea. A total
els and performance, was featured in docking predictions of 98 acceptable models were submitted for T92, of
for the other targets in this category (see Supporting which only 12 were of medium quality, but the models
Information Tables S2 and S3). The models of individual were contributed by a large fraction of the participating
subunits were also of similar accuracy or higher. groups (17 out of 23). On the other hand, the scorer
Excellent performance was obtained for targets T80 performance was very good with 68 acceptable models of
(T0819) and T93 (T0851) with >100 correct models of which almost half (33) were of medium quality. These
which (cid:2)70 were of medium quality, followed by targets models were contributed across most scorer groups (10
T90 (T0843) and T91 (T0847), for which >100 correct
out of 11). CASP participants achieved a particularly
models, comprising (cid:2)40 medium quality ones- were
good performance. Of the 23 models submitted by CASP
submitted. These targets featured subunits sizes of 176–
groups, 17 were of acceptable quality or better, and those
456 residues.
were contributed by six of the seven participating groups.
T80 (T0801) was the sugar aminotransferase WecE
The accuracy of the subunit models was in general lower,
from E.coli K-12, with 376 residues per subunit. Submit- with LGA_S (cid:2)70 and rmsd (cid:2)7 A˚ for the best models,
ted models were evaluated against one interface (1960
and LGA_S values of 50 – 60 for most other models.
A˚2) between the two subunits of the crystal asymmetric
In T94, predicted complexes were assessed only against
unit [Fig. 2(c)]. A total of 27 CAPRI predictor groups the largest interface (1190 A˚2), formed between large
submitted 105 models of acceptable quality or better.
domains of the adjacent subunits, as the second largest
The majority of these (71 models) were of medium qual- interface was much smaller (620 A˚2). In total, 97 accept-
ity. 12 CAPRI groups participated in the scoring experi-
able homodimer models only, were contributed for this
ment and submitted 120 models, of which about half
target: 58 models by CAPRI predictors, 37 by CAPRI
(51) were of medium quality and 14 were acceptable
scorers, and 2 by CASP groups [see Supplementary Table
models. Six CASP participants submitted 11 medium
S2, and Fig. 2(g,h) for a pictorial summary]. The lower
quality models, and two models of acceptable quality.
accuracy of the subunit models for this target (LGA_S
The top ranking CAPRI predictor groups for this target
score (cid:2)58 and rmsd >6 A˚, for the best model) may have
were those of Sali, Guerois, and Eisenstein who submit-
limited the accuracy of the modeled complexes, without
ted 10 medium quality models each. These three groups
however compromising the task of achieving correct
were closely followed by the groups of Seok, Zou, Shen
solutions.
and Lee, each of whom predicted at least five medium
quality models. Each of the three participating servers,
Difficultorproblematichomodimertargets:
HADDOCK, GRAMM-X, and CLUSPRO, submitted at
T68,T72,T77,T79,T86,T88
least one acceptable model. The best performers from
among the scorer groups were those of Zou and Huang This category comprises 6 targets, representing partic-
with 10 medium quality models each, followed by Gray, ular challenges to docking calculations for reason inher-
Kihara and Weng with at least 5 medium quality models, ent to the proteins involved, or targets for which the
and by Fernandez-Recio and Bates with four medium oligomeric state was probably assigned incorrectly at the
quality models. The global landscape of the predictions time of the experiment.
for this target is shown in Figure 2(d). With the exception of T72, targets in this category are
The subunit models for this target were of very high much smaller proteins, than those of the easy dimer tar-
quality, with the best models featuring a LGA_S score of gets (Table I). In three of the targets (T68, T79, T86) the
(cid:2)95 and a backbone rmsd of 1.3 A˚. The quality of the largest interface area between subunits in the crystal is
best models for targets T90 and T91 for which a simi- small (470–860 A˚2) and their oligomeric state assign-
larly high performance was achieved was only somewhat ments were often ambiguous. In the following, we com-
lower, with LGA_S values of 70–88 and backbone rmsd ment on the insights gained from the results obtained
of 2.0–5.0 A˚. for several of these targets.
Interestingly, T91 (T0847), the human Bj-Tsa-9, was No acceptable homodimer models were contributed by
predicted to be a dimer by PISA, but assigned as a CAPRI or CASP groups for targets T68, T77 and T88.
PROTEINS 333

M.F.Lensinketal.
C
O
L
O
R
Figure 3
Targetstructuresandpredictionresultsfordifficultorproblematicdimertargets.T68(T0759),Plectin1and2RepeatsoftheHumanPeriplakin,
PDBcode4Q28.(a)Targetstructureincartoonrepresentation,displaying4subunitsinthecrystal.TheHis-Tagsequence,highlightedinblack,
mediatescontactsatthelargestinterface.(b)Globaldockingpredictionresultsdisplayingonesubunitincartoonrepresentation,withthecenterof
massofthesecondsubunitinthetarget(redsphere),andindockingsolutionssubmittedbyCAPRIpredictors(lightbluespheres),CAPRIscorers
(darkbluespheres),andCASPpredictors(yellowspheres).T77(T0780),conservedhypotheticalprotein(SP_1560),Streptococcuspneumoniae
TIGR4PDBcode4QDY.(c)Targetstructure,highlightingtheassessedinterface(dashedline).(d)Globaldockingpredictionresultsbydifferent
predictorgroups(seelegend(b)fordetail).T88(T0825),syntheticwrapfiveprotein(structureunreleased).(e)Targetstructure.(f)Globaldock-
ingpredictionresultsbydifferentpredictorgroups.T72(T0772),SusDhomolog(BT2259)fromBacteroidesthetaiotaomicronVPI-5482,PDBcode
4Q69.(g)Targetstructure,highlightingthethreeassessedinterfaces.(h)Globaldockingpredictionresultsforthethreeinterfaces,bydifferentpre-
dictorgroups.T79(T0792),OSKAR-N,PDBcode5a49.(i)Targetstructure,highlightingthethreeassessedinterfaces.(j)Globaldockingpredic-
tionresultsforthethreeinterfacesbydifferentpredictorgroups.T86(T0815)Putativepolyketidecyclase(proteinSMa1630)fromSinorhizobium
meliloti,PDBcode4U13.(k)Targetstructure,showingthreeinterfaces.(l)Globaldockingpredictionresultsforthetwointerfacesbydifferentpre-
dictorgroups(theinterfacewiththeyellowmonomerwasnotassessed).
The main problem with T68 (T0759), the plectin 1 and largest interface (860A˚2), but not against the 2 much
2 repeats of the Human Periplakin, was that the crystal smaller interfaces (240 and 160 A˚2).
structure contains an artificial N-terminal peptide repre- Most predictor groups (from both CASP and CAPRI)
senting the His-tag (MGHHHHHHS...) that was used carried out docking calculations without the His-tag,
for protein purification. The N-terminal segments of which they assumed was irrelevant to dimer formation
neighboring subunits, which contain the artificial pep- in-vivo. They were therefore unable to obtain docking
tide, associate to form the largest interface between the solutions that were sufficiently close to the largest inter-
F3 subunits in the crystal (1150 A˚2) [Fig. 3(a)]. Submitted face of the target [Fig. 3(b)]. As well, no acceptable solu-
model were assessed against this interface and the second tions were obtained for second largest interfaces,
334 PROTEINS

PredictionofHomoandHeteroproteinComplexesbyProteinDockingandModeling
indicating that it too was unlikely to represent a stable subunit models for the more truncated subunit were
homodimer. much poorer (rmsd 6.5–10 A˚), and since the helical
The quality of the subunit models was also lower than region of the shorter subunit contributes significantly to
for many other targets (the best model had an LGA_S the dimer interface, whose total area is not very large
score of (cid:2)57), as most groups ignored the His-Tag in ((cid:2)1300 A˚2), no acceptable docking solutions were
building the models as well (see Fig. 6 and CAPRI web- obtained [Fig. 3(e,f)].
site for details). Considering that the His-Tag containing For the other three targets in this category, T72, T79,
peptide contributes significantly to the largest subunit and T86, the homodimer prediction performance
interface, the protein is likely a monomer in absence of remained rather poor, with only very few acceptable
the artificial peptide. This is in fact the authors’ assign- models submitted. The main issue with T79 (T0792), the
ment in the corresponding PDB entry (4Q28), and in OSKAR-N protein, and T86 (T0815), the polyketide
retrospect this target should not have been considered Cyclase from Sinorhizobium meliloti, was likely their very
for the CAPRI docking experiments. small subunit interface (Table I). T79 was predicted by
Different factors contributed to the failure of produc- PISA to be a dimer, but the area of its largest subunit
ing acceptable docking solution for T77 (T0780), the interface is only 680 A˚2. T86, predicted to be dimeric by
conserved hypothetical protein (SP-1560), from Strepto- both PISA and the authors (as stated in the PDB entry,
coccus pneumonia TGR4 [Fig. 3(c,d)]. The protein con- 4U13), has even smaller size subunit interfaces with the
sists of two YbbR-like structural domains (according to largest one burying no >470 A˚2. In both cases these
Pfam) arranged in a crescent-like shape. The domains interfaces are much smaller than the average size
adopt rather twisted b-sheet conformations with exten- required in order to stabilize weak homodimers. 46 It is
sive stretches of coil, and are connected by a single poly- therefore likely that these two proteins are in fact mono-
peptide segment, suggesting that the protein displays an meric at physiological concentrations. Furthermore, T79
appreciable degree of flexibility both within and between and T86 are quite small proteins (80 residues for T79,
the domains. Probably as a consequence of this flexibil- and 100 residues for T86), and it is not uncommon that
ity, the structures of most templates identified by predic- proteins of this size cannot form large enough interfaces
47
tor groups (which approximated only one domain), were unless they are intertwined.
not close enough to that of the target (Supporting Infor- Thisnotwithstanding,afewacceptablehomodimermod-
mation Table S5). As a result, the subunit models were elswerecontributedforallthreeassessedinterfaces(interfa-
generally quite poor, with the best model featuring an ces1,2,3)ofT79(SupportingInformationTableS2).
LGS-A score of only (cid:2)40 (rmsd (cid:2)7 A˚). Although the Among predictor groups, 17 acceptable docking solu-
largest interface of the target is of a respectable size tions (of which five were medium quality models) were
(1600 A˚2) and involves intermolecular contacts between obtained for the largest interface (interface 1). Twelve
one of the domains only, the docking calculations were acceptable solutions, of which one medium quality one,
unable to identify it. The best docking model was incor- were obtained for the second smaller interface (440 A˚2),
rect as it displayed an L-rms (cid:2)19 A˚, and an I-rms (cid:2)10 A˚ and no acceptable quality solutions were obtained for the
(see Supporting Information Table S4). thirdassessedinterface(400A˚2)[seeFig.3(i,j)foranover-
A very different issue plagued the docking prediction view of the prediction results]. Seven CAPRI predictor
of T88 (T0825), the wrap5 protein. The information groups, 1 CASP group and one server (GRAMM-X) con-
given to predictors was that the protein is a synthetic tributed the correct models for interface 1, and seven
construct built from 5 sequence repeats, and is similar to CAPRIgroupssubmittedacceptablemodelsforinterface2.
2YMU (a highly repetitive propeller structure). It was Interestingly scorers did less well than predictors for
furthermore stated that the polypeptide has been mildly interface 1, but better for interface 2, and two scorer
proteolyzed, yielding two slightly different subunits, in groups submitted two acceptable models for interface 3,
which the N-terminus of the first repeat was truncated whereas none were submitted by predictor groups.
to different extent, and that therefore the dimer forms in Overall, the models for the T79 subunit were quite
a non-trivial way. Predictors were given the amino acid accurate, with the best model having and LGA_S score of
sequence of the two alternatively truncated polypeptides. (cid:2)89 and rmsd (cid:2)1.9 A˚.
It turned out that the longer of the two chains, with Not too surprisingly, the dimer prediction perform-
the nearly intact first repeat forms the expected 5-blade ance for T86 was significantly poorer, with only three
b-propeller fold, whereas the chain with the severely acceptable models submitted by CAPRI predictors
truncated first repeat forms only four of the blades, with (Ritchie and Negi) for the largest interface (470 A˚2).
the remainder of the first repeat forming an a-helical Scorers identified five acceptable models for interface 1
segment that contacts the first repeat [Fig. 3(e)]. (Fernandez-Recio and Gray), and two acceptable (or bet-
Both CAPRI and CASP predictor groups were quite ter) models for interface 2 (Seok and Kihara). None of
successful in building very accurate models for the less the 19 models submitted by the seven CASP groups were
truncated subunit (rmsd<0.5 A˚, LGA_S (cid:2)90). But correct [Fig. 3(k,l) for a pictorial summary].
PROTEINS 335

M.F.Lensinketal.
Different problems likely led to the weak prediction Interestingly, acceptable or better models were submit-
performance for Target T72 (T0770), the SusD homolog ted only for the smaller interface (475 A˚2) (Supporting
(BT2259) from Bacteroides Thetaiotaomicron. While the Information Table S2). CAPRI predictors submitted 37
largest subunit interface is of near average size (1120 A˚2), acceptable models, of which 27 were of medium quality,
the interface itself is poorly packed and patchy, an indi- and scorers submitted 27 acceptable models (including
cation that it may not represent a specific association. 21 medium quality ones) [Fig. 4(b)]. Indeed no accepta-
Not too surprisingly, therefore, this led to a poor predic- ble models were submitted for the largest interface (560
tion performance. Overall only three models of accepta- A˚2), which is assigned as the dimer interface in the PDB
ble quality were submitted by CAPRI dockers, namely by entry for this protein.
the HADDOCK and SWARMDOCK servers, and the The failure to model a higher order oligomer for this
Guerois group, each contributing 1 such model. The best target was not due to the quality of the subunit models
of these models (contributed by Guerois) had f(nat) as the latter was quite high (see Fig. 6 and CAPRI web-
(cid:2)29% and L-rms and I-rms values of 8.85 and 3.57 A˚, site), and is probably rooted in the pattern of contacts
respectively. Seven acceptable models were submitted by made by the protein in the crystal, which suggest that
scorers. Bonvin contributed two models, and the groups this target is likely a weak dimer. Considering that all the
of Huang, Grudinin, Gray, Weng and Fernandez-Recio, acceptable docking models involve a different interface
respectively, submitted one model. The best quality mod- than that assigned in the corresponding PDB entry, it is
els had f(nat) (cid:2)18%, and L-rms and I-rms values of furthermore possible that the interface identified in these
(cid:2)7.29 and 4.28 A˚, respectively. No acceptable models solutions is in fact the correct one. But given the very
were submitted by CASP participants. The target struc- small size of either interface, the protein could also be
ture and the distribution of the all the docking solutions monomeric.
A similar situation was encountered with T74
are depicted in Figure 3(g,h).
(T0774), a hypothetical protein from Bacteroides vulga-
The accuracy of the subunit models for T72 was rea-
tus. Here too the target was assigned as a tetramer by
sonable, with the best models having a LGA_S score of
(cid:2)70 (backbone rmsd (cid:2)3.8 A˚). The three successful PISA at the time of the predictions, but is listed as a
monomer by the authors in the PDB entry (4QB7).
CAPRI predictor groups (HADDOCK, SWARMDOCK
Associating the subunits according to the two largest
and Guerois) all had somewhat lower quality subunit
interfaces (520 and 490 A˚2), also produced an open-
models with LGA_S scores in the range of 55 – 67.
ended assembly rather than a closed tetramer, and this
time no acceptable solutions were produced for either
Targetsassignedastetramers:T70,T71,T73,T74,T78
interface, strongly suggesting that the protein is mono-
Five targets were assigned as tetramers at the time of
meric as specified by the authors. It is noteworthy that
the prediction experiment. As described in Assessment
the subunit models for this target were particularly poor
Procedure and Criteria, models for tetramer targets were (LGA_S values (cid:2)40, and rmsd (cid:2)7 A˚), which could also
assessed by systematically comparing all the interfaces in have hampered identifying some of the binding
each model to all the relevant interfaces in the target, interfaces.
and selecting the best-predicted interfaces. Most predic- T71 (T0768), the leucine-rich repeat protein from bac-
tor groups used a two-step approach to build their mod- teroides capillosus, was a difficult case for other reasons.
els. First they derived the model of the most likely dimer, Subunit contacts in the crystal are mediated through
and then docked the dimers to one another. Some three different interfaces, ranging in size from 470 A˚2 to
groups imposed symmetry restraints as part of the dock- 720 A˚2. A closed tetrameric assembly can be built by
ing procedures, or combined this approach with the two- combining interfaces 1 and 3, associating the dimer
step procedure. formed by subunits A and B with the equivalent dimer
In three of the targets (T70, T71, T74) predictors faced of subunits C and D, as shown in Figure 4(c). Interfaces
the problem that all the pair-wise subunit interfaces were 1 and 3 were also those for which some acceptable pre-
quite small (440–720 A˚2), making it difficult to identify dictions were submitted. One acceptable model was con-
stable dimers to initiate the assembly procedure. tributed for the largest interface, by the GRAMM-X, an
T70 (T0765), the modulator protein MzrA from Kleb- automatic server. Eleven acceptable models were submit-
siella Pneumoniae Sub Species, was assigned as a tetramer ted for the third interface (470 A˚2) by 4 CAPRI predictor
at the time of the predictions, but is listed as a dimer groups, and six acceptable models were submitted by
(predicted by PISA and assigned by the authors) in the four CAPRI scorer groups. All the models submitted by
PDB entry (4PWU). Only two of its interfaces in the a single CASP group were wrong. No group succeeded in
F4 crystal bury an area exceeding 400 A˚2 [Fig. 4(a)]. The building the tetramer that comprises the correct models
assembly built by propagating these two interfaces for interfaces 1 and 3 at the same time. Some models
appears to form an extensive layered arrangement across looked promising, but when superimposing equivalent
unit cells in the crystal, rather than a closed tetramer. subunits (in the model vs. the target) the neighboring
336 PROTEINS

PredictionofHomoandHeteroproteinComplexesbyProteinDockingandModeling
C
O
L
O
R
Figure 4
Targetstructuresandpredictionresultsfortetramerictargets.T70(T0765),ModulatorproteinMzrA(KPN_03524)fromKlebsiellapneumoniae
subspecies.(a)Targetstructureincartoonrepresentation,highlightingthetwoassessedinterfaces(dashedlines).(b)Globaldockingprediction
resultsdisplayingonesubunitincartoonrepresentation,withthecenterofmassofthesecondsubunitinthetarget(redspheres),andindocking
solutionssubmittedbyCAPRIpredictors(lightbluespheres),CAPRIscorers(darkbluespheres),andCASPpredictors(yellowspheres).T71
(T0768)Leucine-richrepeatprotein(BACCAP_00569)fromBacteroidescapillosu,PDBcode4QJU.(c)Targetstructureincartoonrepresentation,
highlightingthetworelevantinterfaces(interfaces1and3)(dashedlines).(d)Globaldockingpredictionresultsfortheassessedinterfacesbydif-
ferentpredictorgroups(monomercolorcorrespondingto(c),thatis,theredspheresrepresentthesame,blue,monomer).T73(T0772),Putative
glycosylhydrolase,PDBcode4QHZ.(e)Targetstructureincartoonrepresentation,highlightingthetwoassessedinterfaces(interface1and2)
(dashedlines).(f)Globaldockingpredictionresultsfortheassessedinterfacesbydifferentpredictorgroups.
subunit of the model (the one across the incorrectly pre- None of the predicted tetramer models simultaneously
dicted interface) had its position significantly shifted rel- captured both interfaces, as illustrated in Figure 4(e,f).
ative to that in the target, resulting in an incorrect For T78, no acceptable solutions were submitted by any
structure of the tetrameric assembly. of the participating groups, but the subunit models were
The remaining two targets, T73 (T0772), a putative only marginally more accurate than those of T73.
glycosyl hydrolase from Parabacteroides distaspnos, and The conclusions to be reached from the analysis of
T78 (T0786), a hypothetical protein from Bacillus cereus, these five targets are twofold. One is that the oligomeric
were genuine tetramers assigned as such by both PISA state assignment for higher order assemblies such as tet-
and the authors. Both targets are proteins of similar size ramers is more error prone than that of dimer versus
((cid:2)260 residues) adopting an assembly with classical D monomers. Tetramers often involve smaller interfaces
2
symmetry, which comprises two interfaces, a sizable one between subunits, especially those formed between indi-
(>1000 A˚2) and a smaller one. But the main bottleneck vidual proteins when two dimers associate, and therefore
for both targets was that their larger interface was inter- predictions on the basis of pair-wise crystal contacts such
twined. Available templates did not seem to capture the as those by PISA become unreliable. Independent experi-
intertwined associations, as witnessed from the overall mental evidence is therefore required to validate the exis-
poorer models derived for the individual subunits. For tence of a higher order assembly. The second conclusion
both targets, the best models had an LGA_S score (cid:2)50 to be drawn is that the prediction of higher order assem-
and a backbone rmsd of (cid:2)5–10 A˚. For T73, a total of bly by docking procedures remains a challenge. Accepta-
only nine acceptable models were submitted by the ble models derived for the largest dimer interface are
CAPRI predictor groups of LZERD, Zou and Kihara for probably not accurate enough to enable the identification
the largest interface, and two acceptable models were of stable association modes between two modeled
submitted by the Lee group for the second interface. dimers. This indicates in turn that the propagation of
PROTEINS 337

M.F.Lensinketal.
|     | errors is     | the problem |                  | that currently |            | hampers        | the          | model- |     |     |     |     |     |     |
| --- | ------------- | ----------- | ---------------- | -------------- | ---------- | -------------- | ------------ | ------ | --- | --- | --- | --- | --- | --- |
|     | ing of higher |             | order assemblies |                | from       | the structures |              | of its |     |     |     |     |     |     |
|     | components    |             | in absence       | of             | additional |                | experimental |        |     |     |     |     |     |     |
information.
Heterocomplextargets:T81,T89
|     | T81 (T0797/T0798) |             |                      | and T89   | (T0840/T0841)      |               | were         | the     |     |     |     |     |     |     |
| --- | ----------------- | ----------- | -------------------- | --------- | ------------------ | ------------- | ------------ | ------- | --- | --- | --- | --- | --- | --- |
|     | only two          | bona-fide   | heterocomplex        |           |                    | targets       | in Round     | 30.     |     |     |     |     |     |     |
|     | T81 is            | the complex |                      | between   | the cGMP-dependent |               |              | pro-    |     |     |     |     |     |     |
|     | tein Kinase       | II          | leucine              | zipper    | (44                | residues)     | and          | the     |     |     |     |     |     |     |
|     | Rab11b            | protein     | (198                 | residues) | (PDB               | code          | 4OJK).       | T89 is  |     |     |     |     |     |     |
|     | the complex       |             | between              | the       | much               | larger        | RON receptor |         |     |     |     |     |     |     |
|     | tyrosine          | kinase      | subunit              | (669      | residues)          | and           | the          | macro-  |     |     |     |     |     |     |
|     | phage stimulating |             | protein              | subunit   | (MSP)              | (253          | residues).   |         |     |     |     |     |     |     |
|     | The crystal       |             | structure            | of T81    | features           | two           | Rab11b       | pro-    |     |     |     |     |     |     |
|     | teins binding     |             | on opposite          | sides     | of                 | the centrally |              | located |     |     |     |     |     |     |
|     | leucine           | zipper,     | in a quasi-symmetric |           |                    | arrangement,  |              | which   |     |     |     |     |     | C   |
O
|     | likely represents |     | the | stoichiometry |     | of the | biological | unit |     |     |     |     |     |     |
| --- | ----------------- | --- | --- | ------------- | --- | ------ | ---------- | ---- | --- | --- | --- | --- | --- | --- |
L
| F5  | [Fig. 5(a)]. | A   | total of | 3 interfaces | were | evaluated |     | for this |     |     |     |     |     |     |
| --- | ------------ | --- | -------- | ------------ | ---- | --------- | --- | -------- | --- | --- | --- | --- | --- | --- |
O
|     | targets:    | Interface | 1 (chains |           | C:A, leucine | zipper    | helix       | 1/          |     |     |     |     |     | R   |
| --- | ----------- | --------- | --------- | --------- | ------------ | --------- | ----------- | ----------- | --- | --- | --- | --- | --- | --- |
|     | one copy    | of the    | Rab11b    | protein), |              | Interface | 2 (C:D,     | leu-        |     |     |     |     |     |     |
|     | cine zipper | helix     | 1/helix   | 2),       | interface    | 3         | (equivalent | to Figure 5 |     |     |     |     |     |     |
interface 1). The two Rab11b/zipper helix interfaces were Targetstructuresandpredictionresultsforheterocomplextargets.T81
not exactly identical (780 A˚2 for interface 1 and 630 A˚2 (T0797/T0798),cGMP-dependentProteinKinaseIILeucinZipperand
Rab11bProteinComplex,PDBcode4OJK.(a)Targetstructureincar-
|     | for interface | 2). | The | interface | between | the | helices | of the |     |     |     |     |     |     |
| --- | ------------- | --- | --- | --------- | ------- | --- | ------- | ------ | --- | --- | --- | --- | --- | --- |
toonrepresentation,highlightingtheinterfaceoftheleucinezipper
|     | leucine | zipper | was somewhat |     | larger | (780 | A˚2). | Overall, |     |     |     |     |     |     |
| --- | ------- | ------ | ------------ | --- | ------ | ---- | ----- | -------- | --- | --- | --- | --- | --- | --- |
dimer(2),andthetwoequivalentinterfaces(1,3),betweenthezipper
the interface area of a single copy of the Rab11b protein dimerandthetwoRab11bproteins(dashedlines).(b)Globaldocking
binding to the leucine zipper dimer measures 1070 A˚2. predictionresultsdisplayingoneoftheRab11bsubunitsincartoonrep-
resentation,withthecenterofmassoftheleucinezipperdimerinthe
Consolidating correct predictions for the equivalent target(redsphere),andindockingsolutionssubmittedbyCAPRIpre-
interfaces (Interfaces 1 and 3), the prediction perform- dictors(lightbluespheres),CAPRIscorers(darkbluespheres),and
ance for this complex as a whole was disappointing. CASPpredictors(yellowspheres).T89(T0840/T0841),complexofthe
RONreceptortyrosinekinasesubunitandthemacrophagestimulating
|     | Only 12 | correct | models | were | submitted | by  | the 7 | CAPRI |     |     |     |     |     |     |
| --- | ------- | ------- | ------ | ---- | --------- | --- | ----- | ----- | --- | --- | --- | --- | --- | --- |
proteinsubunit(MSP)(structurenotreleased).(c)Targetstructurein
|     | predictor | groups | of  | Guerois, | Seok, | Huang, | Vajda/Koza- |     |     |     |     |     |     |     |
| --- | --------- | ------ | --- | -------- | ----- | ------ | ----------- | --- | --- | --- | --- | --- | --- | --- |
cartoonrepresentation.(d)Globaldockingpredictionresultsdisplaying
kov, SWARMDOCK, CLUSPRO (a server) and Bates. theRONreceptorkinassubunit,incartoonrepresentations,andthe
Five of those (submitted by Guerois, Seok and Huang) centerofmassoftheMCPproteinsinthetargetandindockingsolu-
tionssubmittedbydifferentpredictorgroups.
|     | were of | medium      | quality. | The        | performance |             | of       | CAPRI |     |     |     |     |     |     |
| --- | ------- | ----------- | -------- | ---------- | ----------- | ----------- | -------- | ----- | --- | --- | --- | --- | --- | --- |
|     | scorers | was better, | with     | 54 correct | models      |             | of which | 16 of |     |     |     |     |     |     |
|     | medium  | quality.    | All      | 11 scorer  | groups      | contributed |          | these |     |     |     |     |     |     |
models, and the best scorer performance was achieved by their performance was much poorer than that of CAPRI
|     | the groups | of  | Bates, | followed | by those | of  | LZERD, | Oliva,       |           |     |            |     |           |     |
| --- | ---------- | --- | ------ | -------- | -------- | --- | ------ | ------------ | --------- | --- | ---------- | --- | --------- | --- |
|     |            |     |        |          |          |     |        | groups. Only | 23 models |     | out of the | 223 | submitted | by  |
Huang, Fernandez-Recio and Seok. The prediction land- CASP groups (10%) were correct, and 6 of these were
|     | scape for | this | target is | shown | in Figure | 5(b). |     |                 |         |     |     |     |     |     |
| --- | --------- | ---- | --------- | ----- | --------- | ----- | --- | --------------- | ------- | --- | --- | --- | --- | --- |
|     |           |      |           |       |           |       |     | medium accuracy | models. |     |     |     |     |     |
T89, the RON receptor kinase subunit complex with The best performance among CAPRI predictor groups
|     | MSP, was    | a simpler      | target, | given | the       | clear, | binary | charac-       |          |           |          |       |            |       |
| --- | ----------- | -------------- | ------- | ----- | --------- | ------ | ------ | ------------- | -------- | --------- | -------- | ----- | ---------- | ----- |
|     |             |                |         |       |           |        |        | was by the    | HADDOCK  | server,   | followed | by    | the groups | of    |
|     | ter of this | heterocomplex. |         | But   | the large | size   | of the | recep-        |          |           |          |       |            |       |
|     |             |                |         |       |           |        |        | Vakser, Seok, | Guerois, | Grudinin, | Lee,     | Huang | and        | Tomii |
tor subunit, and the relatively small interface it formed (see Supporting Information Table S2). A pictorial sum-
|     | with MSP, | represented |     | a challenge | for | the | docking | calcu-  |                |             |     |     |             |     |
| --- | --------- | ----------- | --- | ----------- | --- | --- | ------- | ------- | -------------- | ----------- | --- | --- | ----------- | --- |
|     |           |             |     |             |     |     |         | mary of | the prediction | performance |     | for | this target | is  |
lations. The prediction performance for this complex was provided in Figure 5(c,d).
|     | quite good | overall,       | with | a total      | of 87 | correct | models        | sub- |     |     |     |     |     |     |
| --- | ---------- | -------------- | ---- | ------------ | ----- | ------- | ------------- | ---- | --- | --- | --- | --- | --- | --- |
|     | mitted     | by predictors, |      | representing |       | 41% of  | all submitted |      |     |     |     |     |     |     |
predictor models. Unlike for many other targets of this Resultsacrosstargetsandgroups
|     | round, | scorers | did only | marginally |     | better, | with | 42% of |     |     |     |     |     |     |
| --- | ------ | ------- | -------- | ---------- | --- | ------- | ---- | ------ | --- | --- | --- | --- | --- | --- |
AcrosstargetperformanceofCAPRIdockingpredictions
|     | correct | models. | CASP | groups | were | specifically | invited | to  |     |     |     |     |     |     |
| --- | ------- | ------- | ---- | ------ | ---- | ------------ | ------- | --- | --- | --- | --- | --- | --- | --- |
submit models for this target, and 55 groups did, nearly Results of the docking and scoring predictions for the
ten times more than for other targets in this round. But 25 assessed targets of Round 30, obtained by all groups
|     | 338 PROTEINS |     |     |     |     |     |     |     |     |     |     |     |     |     |
| --- | ------------ | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |

PredictionofHomoandHeteroproteinComplexesbyProteinDockingandModeling
that submitted models for at least one target, are sum- limitations of the docking or modeling procedures. Two
marized in Figure 6 and in the Supporting Information of the targets, T70 and T74, seem to have been errone-
Table S3. For a full account of the results for this Round ously assigned as tetramers at the time of the prediction
the reader is referred to the CAPRI web site (http://www. by PISA, as described above. T70 was assigned as a
ebi.ac.uk/msd-srv/capri/). dimer, and T74 as a monomer, by the respective authors
The results summarized in Figure 6 show clearly that in the PDB entry. In agreement with the authors’ assign-
the prediction performance varies significantly for targets ment, no acceptable solutions were identified for any of
in the four different categories. As expected, the per- the interfaces in T74. Somewhat surprisingly, the quality
formance is significantly better for the 12 dimer targets of the subunits models for this target was particularly
in the “easy” category, than for those in the other cate- poor as well (average LGA_S (cid:2)30).
gories. For 10 of the 12 “easy” targets, at least 30% of In T70, the docking calculations were able to identify
the submitted models per target are of acceptable quality only the smaller of the two interfaces as forming the
or better, and for most of these (eight out of 10), at least dimer interface (Fig. 6), but this interface differs from
20% of the models are of medium quality. The accuracy the one assigned by the authors. This result leaves open
of the subunit models (top panel, Fig. 6) is rather good the possibility that this protein may indeed be a weak
for most of these targets. With the exception of T93, for dimer, in agreement with the author’s assignment, albeit
which the quality of the subunits models spans a wide a different dimer than the one that they propose. Thus
range (LGA_S (cid:2)40–80), the models of the remaining 11 for both of these seemingly erroneously assigned tetra-
targets achieve high LGA_S scores with averages of 80 or meters, the docking calculations actually gave the correct
above. answer, which supports the author’s subsequent assign-
The two less well-predicted targets in this category are ments, which were not made available at the time of the
T92 and T94, probably due to the lower quality of the prediction experiment.
subunit models (average LGA_S<60) (top panel, Fig. 6). For the other three tetrameric targets, T71, T73 and
The docking prediction performance is quite poor for T78, the poor interface prediction performance reflects
the six “difficult or problematic” dimer targets, where a the genuine challenges of modeling higher order oligom-
few acceptable models were submitted for only three of ers. In T71 the small size of the individual interfaces was
the targets (T72, T79, T86), and no acceptable models likely the reason for the paucity of acceptable dimer
were submitted for the remaining three targets. This very models, and those were moreover not accurate enough
poor performance was not rooted in the docking or to enable the correct modeling of the higher order
modeling procedures but rather in the targets themselves. assembly (dimer of dimers). In T73 and T78, the very
In 4 of the targets in this category (T68, T72, T79, T86) few acceptable models for interfaces in the former, and
the oligomeric state (dimer in this case), often predicted the complete failure to model any of the interfaces in the
only by PISA, but sometimes also provided by the latter (Fig. 6), likely stem from the lower accuracy of the
authors, was likely incorrectly assigned. In T68, the His- corresponding subunit models (average LGA_S (cid:2)50–60).
tag used for protein purification and included in the The docking prediction performance was better, but
crystallization forms the observed dimer interface, which not particularly impressive for the two heterocomplex
is therefore most certainly non-native. In T72 the main targets T81 and T89, which represent the type of targets
problem was its very poorly packed and patchy interface, that the CAPRI community commonly deals with. For
suggesting that the dimer might be a crystal artifact, T81 only (cid:2)5% of the submitted models were of accepta-
whereas in T79 and T86, all the pair-wise interfaces in ble quality or better, whereas for T89 the corresponding
the crystal structure were too small for any of them to model fraction was 40%, similar to that achieved for the
represent a stable dimer. easy dimer targets. The poorer performance for T81 can
The only genuinely difficult dimer targets were T77 be readily explained by the fact that this target was in
and T88. For T77, the subunits of this flexible 2-domain fact a hetero tetramer, two copies of the Rab11b protein
protein were rather poorly modeled (average LGA_S 30– binding to opposite sides of a leucine zipper, which had
40), making it difficult to model the “handshake” to be modeled first.
arrangement of the subunits in the dimer [Fig. 3(c,d)]. These results taken together indicate that homology
In T88, the synthetic wrap5 protein, most predictor modeling techniques and docking calculations are able to
groups failed to meet the challenge of correctly modeling predict rather well the structures of biologically relevant
the shorter of the two subunits, in turn leading to incor- homodimers. In addition we see that the prediction per-
rect solutions for the heterodimer. formance for such targets is on average superior than that
As already mentioned, a very poor performance was obtained for heterocomplexes in previous CAPRI rounds,
achieved for the five targets assigned as tetramers at the where on average only about 10–15% of the submitted
time of the predictions. This is illustrated at the level of models are correct for any given target (http://onlinelibrary.
F6 the individual interfaces in these targets (Fig. 6). How- wiley.com/doi/10.1002/9781118889886.ch4/summary), com-
ever, here too the problem was not necessarily rooted in pared to 25% obtained for the majority of the genuine
PROTEINS 339

M.F.Lensinketal.
C
O
L
O
R
Figure 6
PictorialsummaryofthepredictionresultsperassessedinterfaceofthetargetsinCAPRIRound30.Thelowerpaneldepictsthefractionofmodels
ofacceptableandmediumqualityrespectively,submittedbyCAPRIandCASPpredictorgroups,forthe42assessedinterfacesinall25targets
(listedalongthehorizontalaxis).ThedigitfollowingtheCAPRItargetnumberrepresentstheassessedinterface.Thesymmetrytransformationcor-
respondingtotheassessedinterfacesineachtargetarelistedintheSupportingInformationTableS1.Thefractionofcorrectmodelsisshownsepa-
ratelyforthefourmaintargetcategories:Easydimertargets,difficult(orproblematic)dimertargets,tetramerictargets,andheterocomplextargets.
ThemiddlepaneldisplaysthesamedataformodelssubmittedforthesameinterfacesbyCAPRIscorergroups.Thetoppanelshowsboxplotsof
theLGA_Sscorevaluesofthesubunitsinsubmittedmodelsforthetargetslistedalongthehorizontalaxis.TheLGA_SscoreisoneoftheCASP
35
measuresoftheaccuracyofthepredicted3Dstructureofaprotein. ThereddotsrepresenttheLGA_Sscoreofthesubunitstructureofthebest
qualityhomoorheterocomplexmodelsubmittedforeachtarget.ThebestqualitymodelisdefinedastheonewiththelowestI-rms(seeFig.1for
details).
dimer targets in this Round, including both easy and diffi-
AcrosstargetperformanceofCAPRIscoringpredictions
cult homodimers. This result is not surprising, as interfaces
of homodimers are in general larger and more hydrophobic As shown in the middle panel of Figure 6, CAPRI
than those of heterocomplexes, 45 properties which should scorer groups achieved overall a better prediction per-
make them easier to predict. formance than predictor groups. The scoring experiment
Another noteworthy observation is that docking calcu- involves no docking calculations, and only requires sin-
lations can often help to more reliably assign the protein gling out correct solutions from among the ensemble of
oligomeric state, especially in cases where available models uploaded by groups participating in the docking
assignments were ambiguous. Such cases were encoun- predictions. Clearly, such solutions cannot be identified
tered for several of the difficult or problematic targets, if the ensemble of uploaded models contains only incor-
and for targets assigned as tetramers. On the other hand, rect solutions. Therefore no correct scoring solutions
the main challenge in correctly modeling tetramers is to were submitted by scorers for targets where no acceptable
minimize the propagation of errors caused by even small docking solutions were present within the 100 models
inaccuracies in modeling individual interfaces, which can uploaded by predictor groups for given target.
in turn be exacerbated by inaccurate 3D models of the However, for targets where at least a few correct dock-
protein components. ing models were obtained by predictors, scorers were
340 PROTEINS

PredictionofHomoandHeteroproteinComplexesbyProteinDockingandModeling
often able to identify a good fraction of these models, as high quality models is ranked higher, and when two
well as other models that were not identified amongst the groups submitted the same number of high quality mod-
10 best models by the groups that submitted them (Fig. els, the group with more acceptable models is ranked
6). This was particularly apparent for the easy dimer tar- higher.
gets, where scorers often submitted a significantly higher Overall, a total of 11 CAPRI predictor groups submit-
fraction of acceptable-or-better models (>50%) than in ted correct models for at least 10 targets, and medium
the docking experiment, where this fraction rarely quality models for at least seven targets. These groups
exceeded 40%. A similar result was achieved for the heter- submitted models for at least 20 of the targets. Among
ocomplexes, and was particularly impressive for T81, those, the highest-ranking groups in this Round are
where nearly half of the submitted models by scorers were Seok, Huang, and Guerois, with correct models for 15 or
correct, compared to only 5% for the docking predictions. 16 targets, and medium quality models for 12–14 of
The seemingly superior performance of scorers over these targets. These are followed by Zou, Shen and Gru-
dockers has been observed in previous CAPRI assess- dinin (correct models for 11–14 targets, and medium
ments 16,19 where it was attributed in part to the gener- quality models for 10 or 11 of those). The remaining five
ally poor ranking of models by predictors. Their highest- highest ranking groups, Weng, Vakser, Vajda/Kozakov,
ranking models are often not the highest-quality models, Fernandez-Recio and Lee, achieve correct predictions for
and acceptable or better models can often be found 10–15 targets and medium quality predictions for 7–9 of
those. It is noteworthy that two of the three top ranking
lower down the list and amongst the 100 uploaded mod-
predictor groups (Seok and Guerois), and at least one
els. Another reason is the fact that the search space that
other group (Vakser) made heavy use of template-based
scorers have to deal with is orders of magnitude smaller
modeling, an indication that this approach can be quite
(a few thousands of models), than the search space dock-
effective.
ers commonly sample (tens of millions of models). This
The remaining groups listed in Table IV were ranked
significantly increases the odds of singling correct solu-
lower, as they corrected predicted between 1 and 8 tar-
tions in the scoring experiment.
gets only, and produced only a few medium quality
Clearly however, there is more to the scorers’ perform-
models for these targets. However some of these groups
ance than chance alone, particularly in this CAPRI Round,
submitted predictions for a smaller number of targets.
where the main challenge was to model homo-oligomers.
Their performance can therefore not be fairly compared
Some groups that have also implemented docking servers
to that of other groups.
had their server perform the docking predictions com-
Of the 6 CAPRI automatic docking servers ranked in
pletely automatically, but carried out the scoring predic-
Table IV, HADDOCK and CLUSPRO rank highest, fol-
tions in a manual mode, which still tends to be more
lowed by SWARMDOCK, and GRAMM-X.
robust. In addition, a meta-analysis of the uploaded mod-
It is interesting to note that two top ranking CAPRI
els, such as clustering similar docking solutions and select-
servers submitted correct predictions for 16 targets, just
ing and refining solutions from the most populate clusters
as many as the top ranking predictor groups. But the lat-
can also lead to improved performance.
ter groups still produce more medium accuracy models
This notwithstanding, the actual scoring functions
(>10) than the servers (no more than 9). Thus as
used by scorer groups must play a crucial role. But this
already noted in previous CAPRI assessment, some
role is currently difficult to quantify in the context of
CAPRI servers perform nearly on par with more manual
this assessment.
predictions.
Among the CASP predictor and server groups listed in
PerformanceacrossCAPRIandCASPpredictors,
Table IV, the groups of Umeyama and Dunbrack rank
scorersandservers
highest, and both would rank among the best CAPRI
The ranking of CAPRI-CASP11 participants by their predictor groups as their success rate (fraction of correct
prediction performance on the 25 targets of Round 30 is over submitted models) was also high. Of the servers,
T4 summarized in Table IV. The per-target ranking and per- ROSETTASERVER and SEOK_SERVER rank highest,
formance of participants can be found in the Supporting with a performance level similar to SWARMDOCK.
Information Tables S2 and S3. Thirty-nine CASP groups submitting models for 1–5 tar-
The ranking in Table IV considers only the best quality gets, none of which were correct, are not explicitly listed
model submitted by each group for every target. The in the Table.
ranking in the Supporting Information Table S2 takes Lastly, judging also by the best model submitted for
into account both the total number of acceptable models, each target, CAPRI scorers outperform CAPRI predic-
and the number of higher quality models (medium qual- tors, as already mentioned when analyzing the perform-
ity ones for this Round, as detailed in the section on ance across targets. Highly ranking scorer groups
assessment criteria). When two groups submitted the submitted on average correct models for 1–2 more tar-
same number of acceptable models, the one with more gets than CAPRI predictors, and the number of medium
PROTEINS 341

M.F.Lensinketal.
TableIV quality models that groups submit for these targets is
ParticipantrankingbyTargetperformanceParticipant also somewhat higher.
|     |     |     |     | Of the | 13 scorer | groups | that submitted |     | an accurate |
| --- | --- | --- | --- | ------ | --------- | ------ | -------------- | --- | ----------- |
Participated
targets Performance model for at least one target, 11 have correctly predicted
|     |     |     |     | at least | 10 targets | and submitted | medium | quality | models |
| --- | --- | --- | --- | -------- | ---------- | ------------- | ------ | ------- | ------ |
CAPRIPredictorRanking
| Seok |     | 25  | 15/14** | for seven | of those. |     |     |     |     |
| ---- | --- | --- | ------- | --------- | --------- | --- | --- | --- | --- |
Huang 25 16/13** The best performing groups are those of Bonvin,
| Guerois |     | 25  | 16/12** |               |        |                |           |          |            |
| ------- | --- | --- | ------- | ------------- | ------ | -------------- | --------- | -------- | ---------- |
|         |     |     |         | Bates, Huang, | Seok,  | Zou and        | Kihara,   | followed | closely by |
| Zou     |     | 25  | 14/11** |               |        |                |           |          |            |
|         |     |     |         | four other    | groups | that correctly | predicted | at least | 13 tar-    |
| Shen    |     | 25  | 13/11** |               |        |                |           |          |            |
Grudinin 24 11/10** gets, and produced medium quality models for at least
| Weng            |     | 25  | 13/9** |             |        |      |     |     |     |
| --------------- | --- | --- | ------ | ----------- | ------ | ---- | --- | --- | --- |
|                 |     |     |        | 10 of these | (Table | IV). |     |     |     |
| Vakser          |     | 25  | 11/9** |             |        |      |     |     |     |
| Vajda/Kozakov   |     | 24  | 15/8** |             |        |      |     |     |     |
| Fernandez-Recio |     | 25  | 11/8** |             |        |      |     |     |     |
Factorsinfluencingtheprediction
| Lee |     | 20  | 10/7** |     |     |     |     |     |     |
| --- | --- | --- | ------ | --- | --- | --- | --- | --- | --- |
performance
| Tomii      |     | 20  | 8/6** |        |                |       |          |              |         |
| ---------- | --- | --- | ----- | ------ | -------------- | ----- | -------- | ------------ | ------- |
| Sali       |     | 12  | 6/4** |        |                |       |          |              |         |
|            |     |     |       | Unlike | in previous    | CAPRI | rounds,  | Round        | 30 com- |
| Negi       |     | 25  | 7/3** |        |                |       |          |              |         |
|            |     |     |       | prised | solely targets | where | both the | 3D structure | of the  |
| Eisenstein |     | 6   | 3**   |        |                |       |          |              |         |
Bates 25 7/2** protein subunits and their association modes had to be
| Kihara |     | 23  | 7/2** |          |          |            |             |     |             |
| ------ | --- | --- | ----- | -------- | -------- | ---------- | ----------- | --- | ----------- |
|        |     |     |       | modeled. | Deriving | the atomic | coordinates |     | of the pre- |
| Zhou   |     | 25  | 4/2** |          |          |            |             |     |             |
Tovchigrechko 12 3/1** dicted homo-oligomers therefore involved a number of
Ritchie 8 2/1** steps each requiring the use of specialized software and
| Fernandez-Fuentes |     | 14  | 1   |              |           |             |           |            |          |
| ----------------- | --- | --- | --- | ------------ | --------- | ----------- | --------- | ---------- | -------- |
|                   |     |     |     | making       | strategic | choices as  | to how it | should be  | applied. |
| Xiao              |     | 11  | 1   |              |           |             |           |            |          |
|                   |     |     |     | As mentioned |           | in Synopsis | of the    | Prediction | Methods, |
| Gong              |     | 8   | 0   |              |           |             |           |            |          |
DelCarpio 3 0 the approaches for modeling the subunit structures and
| Wade |     | 2   | 0   | generating | the oligomer |     | models vary | widely | amongst |
| ---- | --- | --- | --- | ---------- | ------------ | --- | ----------- | ------ | ------- |
Haliloglu 1 0 predictor groups, and across targets. It is therefore diffi-
CAPRISERVERRanking
HADDOCK 25 16/9** cult to reliably pinpoint specific factors that contributed
CLUSPRO 25 16/8** or hampered successful predictions. Nonetheless some
| SWARMDOCK |     | 25  | 11/4** |           |            |               |          |        |            |
| --------- | --- | --- | ------ | --------- | ---------- | ------------- | -------- | ------ | ---------- |
|           |     |     |        | general   | trends can | be outlined.  | Even     | though | Round 30   |
| GRAMM-X   |     | 22  | 6/1**  |           |            |               |          |        |            |
|           |     |     |        | comprised | only       | targets whose | subunits | could  | be readily |
| LZERD     |     | 25  | 3      |           |            |               |          |        |            |
DOCK/PIERR 2 1 modeled using templates from the PDB, the subunit
CAPRIScorerRanking modeling strategy had an important influence on the
| Bonvin |     | 25  | 18/14** |                |         |        |      |              |           |
| ------ | --- | --- | ------- | -------------- | ------- | ------ | ---- | ------------ | --------- |
|        |     |     |         | final oligomer | models. | Groups | that | used several | different |
| Bates  |     | 24  | 17/13** |                |         |        |      |              |           |
Huang,Seok 25 16/13** subunit models for the same target increased their
| Zou,Kihara      |     | 25  | 15/12** |                |             |          |               |           |            |
| --------------- | --- | --- | ------- | -------------- | ----------- | -------- | ------------- | --------- | ---------- |
|                 |     |     |         | chance         | of deriving | at least | an acceptable | oligomer  | model.     |
| Fernandez-Recio |     | 25  | 14/12** |                |             |          |               |           |            |
|                 |     |     |         | Such different | models      | were     | obtained      | either by | using dif- |
| Weng            |     | 25  | 16/11** |                |             |          |               |           |            |
Oliva 22 14/11** ferent templates (some groups used as many as five tem-
Grudinin 25 13/10** plates for the same target), or by starting from the same
Gray 17 10/7** template and modifying it by optimizing loop conforma-
| LZERD |     | 25  | 6**   |           |            |              |              |      |             |
| ----- | --- | --- | ----- | --------- | ---------- | ------------ | ------------ | ---- | ----------- |
|       |     |     |       | tions and | subjecting | it to energy | refinements. |      | These opti- |
| Lee   |     | 5   | 3/2** |           |            |              |              |      |             |
|       |     |     |       | mizations | seemed     | particularly | effective    | when | carried out |
| Sali  |     | 1   | 0     |           |            |              |              |      |             |
CASPPredictorandServerRanking
|               |     |     |        | in the context | of             | the oligomers | representing    | the | highest- |
| ------------- | --- | --- | ------ | -------------- | -------------- | ------------- | --------------- | --- | -------- |
| Umeyama       |     | 19  | 13/8** |                |                |               |                 |     |          |
|               |     |     |        | ranking        | template-based | or            | docking models. |     |          |
| ROSETTASERVER |     | 13  | 9/8**  |                |                |               |                 |     |          |
Dunbrack 12 11/6** As already mentioned, information on oligomeric tem-
SEOK_SERVER 22 7/5** plates in the PDB was another important element con-
| Luethy   |     | 8   | 5/4** |             |            |          |            |              |             |
| -------- | --- | --- | ----- | ----------- | ---------- | -------- | ---------- | ------------ | ----------- |
|          |     |     |       | tributing   | to improve | the      | prediction | performance. | This        |
| Nakamura |     | 12  | 7/3** |             |            |          |            |              |             |
|          |     |     |       | information | was        | the main | ingredient | for two      | of the best |
| Baker    |     | 8   | 3**   |             |            |          |            |              |             |
Wallner 2 1** performing groups that heavily relied on template-based
Skwark,Lee,RAPTOR-X_Wang, 1–4 1 docking. Other groups that performed well used mainly
NNS_Lee
|                         |     |     |     | ab-initio | docking          | methods | of various | origins,    | but either |
| ----------------------- | --- | --- | --- | --------- | ---------------- | ------- | ---------- | ----------- | ---------- |
| 39participantsnotlisted |     | 1–5 | 0   |           |                  |         |            |             |            |
|                         |     |     |     | guided    | the calculations | or      | filtered   | the results | based on   |
Foreachtargetonlythebestqualitysolutioniscounted;intotal25targetswere structural information from homologous oligomers.
| assessed. Column | 2 indicates the number | of targets for which | predictions were |       |           |           |      |              |           |
| ---------------- | ---------------------- | -------------------- | ---------------- | ----- | --------- | --------- | ---- | ------------ | --------- |
|                  |                        |                      |                  | Other | important | elements, | such | as selecting | represen- |
submitted.InColumn3,thenumberswithoutstarsindicatemodelsofacceptable
|     |     |     |     | tative members | of  | clusters | of docking | solutions, | and the |
| --- | --- | --- | --- | -------------- | --- | -------- | ---------- | ---------- | ------- |
qualityorbetter,andthenumberswith“**”indicatethenumberofthosemodels
thatwereofmediumquality.
|     |     |     |     | final scoring | functions | used | to rank | models | and select |
| --- | --- | --- | --- | ------------- | --------- | ---- | ------- | ------ | ---------- |
342 PROTEINS

PredictionofHomoandHeteroproteinComplexesbyProteinDockingandModeling
|     |     |     |     |     |     |     |     | complexes      | for       | the 25        | targets     | in this    | Round      | are        | plotted in |     |
| --- | --- | --- | --- | --- | --- | --- | --- | -------------- | --------- | ------------- | ----------- | ---------- | ---------- | ---------- | ---------- | --- |
|     |     |     |     |     |     |     |     | Figure         | 7 as      | a function    | of          | the I-rms  | value.     | The        | LGA_S      |     |
|     |     |     |     |     |     |     |     | measure        | was       | used because  | it          | does not   | depend     | on         | the res-   |     |
|     |     |     |     |     |     |     |     | idue numbering |           | along         | the chain,  | which      | may        | vary       | at least   |     |
|     |     |     |     |     |     |     |     | in a fraction  |           | of the models |             | submitted  | by         | CAPRI      | partici-   |     |
|     |     |     |     |     |     |     |     | pants.         | The I-rms | measure       | was         | used       | as it      | represents | best       |     |
|     |     |     |     |     |     |     |     | the accuracy   |           | level of the  | predicted   |            | interface. |            |            |     |
|     |     |     |     |     |     |     |     | Each           | point     | in Figure     | 7           | represents |            | one        | submitted  |     |
|     |     |     |     |     |     |     |     | model,         | and       | points are    | colored     | according  |            | to the     | quality of |     |
|     |     |     |     |     |     |     |     | the predicted  |           | complex       | (incorrect, | acceptable |            | and        | medium     |     |
|     |     |     |     |     |     |     |     | quality).      | The       | plot clearly  | shows       | that       | medium     | quality    | pre-       |     |
| C   |     |     |     |     |     |     |     | dicted         | complexes | (I-rms        | values      | between    | 1          | and        | 3 A˚) tend |     |
O
| L   |          |     |     |     |     |     |     | to be      | associated    | with    | high    | accuracy |        | subunit   | models   |     |
| --- | -------- | --- | --- | --- | --- | --- | --- | ---------- | ------------- | ------- | ------- | -------- | ------ | --------- | -------- | --- |
| O   |          |     |     |     |     |     |     | (LGA_S     | values        | >80).   | We also | see      | that   | predicted | com-     |     |
| R   |          |     |     |     |     |     |     |            |               |         |         |          |        |           | A˚)      |     |
|     |          |     |     |     |     |     |     | plexes     | of acceptable | quality |         | (I-rms   | values | of 2–4    | are      |     |
|     |          |     |     |     |     |     |     | associated | with          | subunit | models  | that     | span   | a wide    | range in |     |
|     | Figure 7 |     |     |     |     |     |     |            |               |         |         |          |        |           |          |     |
|     |          |     |     |     |     |     |     | accuracy   | levels        | (LGA_S  | between | 30       | and    | 90). This | range    |     |
Subunitmodelaccuracyandthequalityofpredictedcomplexesin is comparable to the subunit accuracy range associated
CAPRIRound30.TheCASPLGA_Sscoresofsubunitmodelsinthe with incorrect models of complexes (I-rms >4 A˚; see
predictedcomplexesforthe25targetsinthisRound(verticalaxis)are
|     |     |     |     |     |     |     |     | Table | III for | details on | how | I-rms | contributes |     | to rank |     |
| --- | --- | --- | --- | --- | --- | --- | --- | ----- | ------- | ---------- | --- | ----- | ----------- | --- | ------- | --- |
plottedasafunctionoftheI-rmsvalues(horizontalaxis).Eachpoint
|     |     |     |     |     |     |     |     | CAPRI | models). | Identical | trends | are | observed | when | plot- |     |
| --- | --- | --- | --- | --- | --- | --- | --- | ----- | -------- | --------- | ------ | --- | -------- | ---- | ----- | --- |
inthisFigurerepresentsonesubmittedmodel,andpointsarecolored
accordingtothequalityofthepredictedcomplex,respectively,incorrect
|     |     |     |     |     |     |     |     | ting the | GDT-TS | scores | as a | function | of  | the I-rms | values |     |
| --- | --- | --- | --- | --- | --- | --- | --- | -------- | ------ | ------ | ---- | -------- | --- | --------- | ------ | --- |
(yellow),acceptable(blue)andmedium(green)quality(seeTableIand
|     |     |     |     |     |     |     |     | for the | fraction | of the | models | with | correct | residues | num- |     |
| --- | --- | --- | --- | --- | --- | --- | --- | ------- | -------- | ------ | ------ | ---- | ------- | -------- | ---- | --- |
thetextfordetails).
|     |             |               |             |        |          |                  |            | bering             | (Supporting | Information      |                | Fig.     | S2).    |                  |           |     |
| --- | ----------- | ------------- | ----------- | ------ | -------- | ---------------- | ---------- | ------------------ | ----------- | ---------------- | -------------- | -------- | ------- | ---------------- | --------- | --- |
|     |             |               |             |        |          |                  |            | That               | both        | accurate         | and inaccurate |          | subunit | models           | are       |     |
|     |             |               |             |        |          |                  |            | associated         |             | with incorrectly |                | modeled  |         | complexes        | is        |     |
|     | those to    | be submitted, | also        | played | a role   | as already       | men-       |                    |             |                  |                |          |         |                  |           |     |
|     |             |               |             |        |          |                  |            | expected.          | Inaccurate  | subunit          |                | models   | may     | indeed           | prevent   |     |
|     | tioned here | and           | in previous | CAPRI  | reports. | 19               |            |                    |             |                  |                |          |         |                  |           |     |
|     |             |               |             |        |          |                  |            | the identification |             | of the           | correct        | binding  | mode,   |                  | and dock- |     |
|     | In the      | following     | we examine  |        | in more  | detail           | the impact |                    |             |                  |                |          |         |                  |           |     |
|     |             |               |             |        |          |                  |            | ing calculations   |             | may              | fail to        | identify | the     | correct          | binding   |     |
|     | of two      | important     | elements    | of     | this     | joint CASP-CAPRI |            |                    |             |                  |                |          |         |                  |           |     |
|     |             |               |             |        |          |                  |            | mode               | even        | when the         | subunit        | models   |         | are sufficiently |           |     |
experiment. We evaluate the influence of the accuracy of accurate. It is however noteworthy that complexes classi-
|     | individual | subunits | models | on  | the oligomer |     | prediction |         |           |        |       |          |     |                 |     |     |
| --- | ---------- | -------- | ------ | --- | ------------ | --- | ---------- | ------- | --------- | ------ | ----- | -------- | --- | --------------- | --- | --- |
|     |            |          |        |     |              |     |            | fied as | incorrect | by the | CAPRI | criteria | do  | not necessarily |     |     |
performance, and estimate the extent to which proce- represent prediction noise, as a recent analysis has shown
|     | dures that | rely        | on docking      | methodology |          | and       | those that |               |     |                 |      |             |             |     |            |     |
| --- | ---------- | ----------- | --------------- | ----------- | -------- | --------- | ---------- | ------------- | --- | --------------- | ---- | ----------- | ----------- | --- | ---------- | --- |
|     |            |             |                 |             |          |           |            | that residues |     | that contribute |      | to the      | interaction |     | interfaces |     |
|     | employ     | specialized | template-based  |             | modeling |           | confer an  |               |     |                 |      |             |             |     |            |     |
|     |            |             |                 |             |          |           |            | are correctly |     | predicted       | in a | significant | fraction    |     | of these   |     |
|     | advantage  | over        | straightforward |             | homology | modeling. |            |               | 48  |                 |      |             |             |     |            |     |
complexes.
|     |     |     |     |     |     |     |     | Somewhat |     | less expected | is  | the observation |     | (Fig. | 7) that | F7  |
| --- | --- | --- | --- | --- | --- | --- | --- | -------- | --- | ------------- | --- | --------------- | --- | ----- | ------- | --- |
Influenceofsubunitmodelaccuracy in a significant number of cases, acceptable and to a
The subunit models used to derive the models of the smaller extent also medium quality docking solutions
|     |           |      |           |        |          |         |       | can be | identified | even | with | lower accuracy |     | models | of the |     |
| --- | --------- | ---- | --------- | ------ | -------- | ------- | ----- | ------ | ---------- | ---- | ---- | -------------- | --- | ------ | ------ | --- |
|     | oligomers | were | generated | either | by CAPRI | groups, | those |        |            |      |      |                |     |        |        |     |
with more homology modeling expertise, or borrowed individual subunits. This is an encouraging observation,
|     |              |     |            |           |     |         |          | as it suggests |     | that docking | calculations |     | can | lead | to useful |     |
| --- | ------------ | --- | ---------- | --------- | --- | ------- | -------- | -------------- | --- | ------------ | ------------ | --- | --- | ---- | --------- | --- |
|     | from amongst |     | the models | submitted |     | by CASP | servers, |                |     |              |              |     |     |      |           |     |
which were made available to CAPRI groups in time for solutions with protein models built by homology, and
|     |                  |           |               |              |          |            |              | that these   | models      | need           | not always   | be         | of the | highest     | accu-     |     |
| --- | ---------------- | --------- | ------------- | ------------ | -------- | ---------- | ------------ | ------------ | ----------- | -------------- | ------------ | ---------- | ------ | ----------- | --------- | --- |
|     | each docking     |           | experiment.   | The          | subunit  | structures | in           |              |             |                |              |            |        |             |           |     |
|     |                  |           |               |              |          |            |              | racy.        | What        | probably       | matters      | more       | for    | the success | of        |     |
|     | models submitted |           | by CAPRI      | and          | CASP     | groups     | for all 25   |              |             |                |              |            |        |             |           |     |
|     |                  |           |               |              |          |            |              | docking      | predictions | is             | the accuracy |            | with   | which       | the bind- |     |
|     | targets of       | Round     | 30 were       | assessed     |          | using      | the standard |              |             |                |              |            |        |             |           |     |
|     |                  |           |               |              |          |            |              | ing regions  | of          | the individual |              | components |        | of the      | complex   |     |
|     | CASP GDT_TS      |           | and LGA_S     | scores,      | as       | well       | as the back- |              |             |                |              |            |        |             |           |     |
|     |                  |           |               |              |          |            |              | are modeled, |             | rather than    | the          | accuracy   | of     | the         | 3D model  |     |
|     | bone rmsd        | of        | the submitted |              | model    | versus     | the target   |              |             |                |              |            |        |             |           |     |
|     |                  |           |               |              |          |            |              | considered   | in          | its entirety.  |              |            |        |             |           |     |
|     | structures.      | The       | values        | of these     | measures |            | obtained for |              |             |                |              |            |        |             |           |     |
|     | models           | submitted | by all        | participants |          | in Round   | 30 for       |              |             |                |              |            |        |             |           |     |
Round30predictionsversusstandardhomologymodeling
|     | each target | can | be found | at the | CAPRI | website | together |     |     |     |     |     |     |     |     |     |
| --- | ----------- | --- | -------- | ------ | ----- | ------- | -------- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
with the assessment results for this Round. To estimate the extent to which docking methods or
To gauge the relation between the accuracy of subunit template-based modeling procedures confer an advantage
models and the docking prediction performance, the over straightforward homology modeling, the accuracy of
LGA-S scores of subunit models in the predicted the submitted oligomer models for each target was
|     |     |     |     |     |     |     |     |     |     |     |     |     |     | PROTEINS | 343 |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | -------- | --- | --- |

M.F.Lensinketal.
TableV
Bestavailabletemplatesdetectedbasedonsequence(“Sequence”),experimentalmonomerstructure(“Monomer”),andexperimentaloligomer
structure(“Oligomer”)Target
BesttemplateTM-score(detectedtemplate)
Targetreleased Databasereleased Sequence Monomer Oligomer
T68 May01,2014 April24,2014 0.348(3njd) 0.370(3fse) 0.370(3fse)
T69 May05,2014 April24,2014 0.852(1qlw) 0.852(1qlw) 0.852(1qlw)
T70 May06,2014 April24,2014 0.639(2f06) 0.644(3c1m) 0.652(3tvi)
T71 May07,2014 April24,2014 0.509(2id5) 0.618(3jur) 0.618(3jur)
T72 May08,2014 April24,2014 0.510(3otn) 0.510(3otn) 0.510(3otn)
T73 May09,2014 April24,2014 a 0.554(1hql) 0.554(1hql)
T74 May12,2014 April24,2014 0.340(4jrf) 0.340(4jrf) 0.340(4jrf)
T75 May13,2014 April24,2014 0.880(3rjt) 0.880(3rjt) 0.880(3rjt)
T77 May15,2014 April24,2014 0.393(2xwx) 0.375(4iib) 0.375(4iib)
T78 May20,2014 May17,2014 0.315(3c6c) 0.370(1o0s) 0.403(2f3o)
T79 May23,2014 May17,2014 0.440(2bnl) 0.469(2xig) 0.471(2w57)
T80 June02,2014 May17,2014 0.938(1mdo) 0.943(2fnu) 0.943(2fnu)
T82 June04,2014 May17,2014 0.846(4dn2) 0.846(4dn2) 0.846(4dn2)
T84 June09,2014 May17,2014 0.939(2btm) 0.941(1b9b) 0.941(1b9b)
T85 June10,2014 May17,2014 0.889(3ggo) 0.889(3ggo) 0.889(3ggo)
T86 June11,2014 May17,2014 0.459(4h3u) 0.467(3gzr) 0.470(3hk4)
T87 June13,2014 May17,2014 0.922(3get) 0.922(3get) 0.922(3get)
T90 July03,2014 June06,2014 0.921(4qgr) 0.927(2oga) 0.927(2oga)
T91 July08,2014 June06,2014 0.750(4gel) 0.750(4gel) 0.808(3hsi)
T92 July09,2014 June06,2014 0.785(1tu7) 0.837(3h1n) 0.837(3h1n)
T93 July10,2014 June06,2014 0.896(4a7p) 0.896(4a7p) 0.896(4a7p)
T94 July11,2014 June06,2014 0.655(3gff) 0.655(3gff) 0.655(3gff)
TM-scoreofthetemplatesthathavethehighestTM-scoreamongtop10selectedtemplatesforeachtargetandthePDBIDsofthetemplatesarelisted.
aNoproteinwiththedesiredoligomerstatewasfoundamongthetop100HHsearchentries.
compared to the accuracy of the models build using the ((cid:3)0.7), are the easy targets, whereas difficult targets are
bestoligomertemplatesforthattargetavailableinthePDB those with poorer templates (lower TM-scores). Many of
atthetimeoftheprediction.Onlydimertargets(andtem- the best templates from all three categories were also
plates) were considered, given the uncertaintyof the oligo- detected and used by predictor groups (see Supporting
mericstateassignmentsforsomeofthetetramerictargets. Information Table S5), even though these groups only
Three categories of the best dimeric templates were had sequence information to identify them during the
considered (see Assessment Procedures and Criteria): prediction round.
templates identified on the basis of sequence alignments The accuracy levels of the models built using the three
alone, templates identified by structurally aligning the categories of best templates for each target and the best
target and template monomers, and templates identified models from each of the participating CAPRI predictor
by structurally aligning the target and template oligom- groups submitted for the same target are plotted in Fig-
ers. Only the sequence-based template selection reconsti- ure 8. The model accuracy is measured by the I-rms F8
tutes the task performed by predictors, to whom only value, representing the accuracy level of the predicted
the target sequence was disclosed at the time of the pre- interface in the complex. Each entry in the Figure repre-
diction. The resulting templates thus represent the best sents one model, and for each template category (based
templates available to predictors during the prediction on sequence alignments, on structural alignment of the
Round. Obviously, the structurally most similar tem- monomers and dimers, respectively), up to 10 best mod-
plates could not be identified by predictors, but are con- els are shown per target and colored according to the
sidered here in order to evaluate the advantage, if any, template category.
conferred by such templates over those identified on the Inspection of Figure 8 indicates that models submitted
basis of sequence alignments. by CAPRI predictor groups, a vast majority of which
T5 Table V lists the best templates from each category employed docking methods as part of their protocol,
identified for all dimeric targets of Round 30 and the tend to be of higher accuracy. For most of the easy tar-
corresponding template-target TM-scores. These tem- gets, the 10 models submitted by CAPRI groups more
plates represent those with the highest TM-score among consistently display lower I-rms values then the models
the best 10 templates from each category detected for a built from the best templates. This is the case not only
given targets. Not too surprisingly targets with more for models derived from the sequence-based templates
similar templates, those featuring high TM-scores but also for the most structurally similar templates of
344 PROTEINS

PredictionofHomoandHeteroproteinComplexesbyProteinDockingandModeling
C
O
L
O
R
Figure 8
AccuracyofRound30homodimermodelspredictedbyproteindockingmethodsandtemplate-basedmodelingversusmodelsderivedbystandard
homologymodeling.TheI-rmsvalues,representingtheaccuracylevelofthepredictedinterface,areplotted(verticalaxis)fordifferentmodelsfor
eachtarget(listedonthehorizontalaxisusingtheCAPRItargetidentification).Eachpointrepresentsonemodel.Thebestmodelssubmittedby
individualCAPRIpredictorgroupsarerepresentedbygreentriangles.Theremainingmodelsarethosebuiltinthisstudybystandardhomology
42
modelingtechniques onthebasisofhomodimertemplatesfromthePDB.Upto10bestmodelsareshownpertargetandtemplatecategory(see
text).Modelsbasedontemplatesidentifiedusingsequenceinformation(blacktriangles),modelsbasedstructuralalignmentsofindividualmono-
mers(redlozenges),andthosebasedonstructuralalignmentsoftheentiredimers(bluetriangles).Thetargets(onlydimers)aresubdividedinto
easyanddifficulttargets(seetext).DashedhorizontallinesrepresentI-rmsvaluesdelimitingmodelsofhigh,medium,acceptableandlower(incor-
rect)qualitybyCAPRIcriteria.
the monomer or dimer categories. Considering only the T85, the highest accuracy models were predicted by the
best models for each targets the performance results are group of Seok, who employed specialized template-based
| mode balanced. |     | For seven | out | of the | 12 easy | targets the |          |            |           |     |         |          |     |
| -------------- | --- | --------- | --- | ------ | ------- | ----------- | -------- | ---------- | --------- | --- | ------- | -------- | --- |
|                |     |           |     |        |         |             | modeling | techniques | augmented |     | by loop | modeling | and |
best models overall were submitted by CAPRI partici- refinement. But the accuracy of these models was not
| pants, whereas | for | the | remaining | five | targets | the most |                 |     |         |          |         |         |     |
| -------------- | --- | --- | --------- | ---- | ------- | -------- | --------------- | --- | ------- | -------- | ------- | ------- | --- |
|                |     |     |           |      |         |          | vastly superior | to  | that of | the best | docking | models. |     |
accurate models were those derived from the structurally Lastly,nottoosurprisingly,oligomermodelsbuildusing
most similar template. Overall however, acceptable or thesequence-basedbesttemplatesweregenerallyofinferior
| medium quality |     | models | were | obtained | with | all the |     |     |     |     |     |     |     |
| -------------- | --- | ------ | ---- | -------- | ---- | ------- | --- | --- | --- | --- | --- | --- | --- |
accuracythanmodelsbuiltfromtemplatesofthetwoother
approaches and for nearly all the easy targets. categories. Interestingly, models derived from the most
| On the | other | hand it | is remarkable |     | that | for three of |              |         |       |           |      |               |     |
| ------ | ----- | ------- | ------------- | --- | ---- | ------------ | ------------ | ------- | ----- | --------- | ---- | ------------- | --- |
|        |       |         |               |     |      |              | structurally | similar | dimer | templates | were | not generally |     |
the difficult targets (T72, T79, and T86), the docking more accurate than those derived from the structurally
| procedures | were | able to | produce | acceptable | models, | with |     |     |     |     |     |     |     |
| ---------- | ---- | ------- | ------- | ---------- | ------- | ---- | --- | --- | --- | --- | --- | --- | --- |
mostsimilarmonomers.Thismaystemfromdifferencesin
| one medium     | quality | model      | for        | T79, | whereas | all the |                |                    |       |      |          |                 |          |
| -------------- | ------- | ---------- | ---------- | ---- | ------- | ------- | -------------- | ------------------ | ----- | ---- | -------- | --------------- | -------- |
|                |         |            |            |      |         |         | the structural | alignmentsthatwere |       |      | used     | to detect       | the tem- |
| template-based | models  | were       | incorrect. |      |         |         |                |                    |       |      |          |                 |          |
|                |         |            |            |      |         |         | plates, which  | in turn            | could | have | affected | the performance |          |
| Overall        | these   | results do | confirm    | that | protein | docking |                |                    |       |      |          |                 |          |
ofthehomologymodelingprocedure(MODELLER).
| procedures     | represent  | an       | added      | value over | straightforward |           |            |     |     |         |     |     |     |
| -------------- | ---------- | -------- | ---------- | ---------- | --------------- | --------- | ---------- | --- | --- | ------- | --- | --- | --- |
| template-based | modeling.  |          | One        | must       | recall however, | that      |            |     |     |         |     |     |     |
| docking        | was often  | combined |            | with       | template-based  |           |            |     |     |         |     |     |     |
|                |            |          |            |            |                 |           | CONCLUDING |     |     | REMARKS |     |     |     |
| restraints     | and hence, | can      | in general |            | not be          | qualified | as         |     |     |         |     |     |     |
ab-initio docking in the context of this experiment. It is CAPRI Round 30, for which results were assessed here,
also important to note that for two targets, T82 and was the first CASP-CAPRI experiment, which brought
|     |     |     |     |     |     |     |     |     |     |     |     | PROTEINS | 345 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | -------- | --- |

M.F.Lensinketal.
together the community of groups developing methods assignment was ambiguous or inaccurate. Such ambigu-
for protein structure prediction and model refinement, ous or inaccurate oligomeric state assignments repre-
with groups developing methods for predicting the 3D sented a confounding factor for the docking prediction
structure of protein assemblies. The 25 targets of this in this round. The problem arose mainly from the fact
round represented a subset of the targets submitted for that the authors’ assignments, usually based on inde-
the CASP11 prediction season of the summer of 2014. In pendent experiment evidence, were not available to pre-
line with the main focus of CASP, the majority of these dictors at the time of the prediction experiment. Instead,
targets were single protein chains, forming mostly homo- predictors were provided with tentative assignments,
dimers, and a few homotetramers. Only two of the tar- inferred on the basis of computational analysis of the
gets were heterodimers, similar to the staple targets in crystal contacts. Quite encouragingly, for most targets
previous CAPRI rounds. Unlike in most previous CAPRI with ambiguous assignment, or for which the tentative
rounds both subunit structures and their association assignments were later overruled by the authors upon
modes had to be modeled for all the targets. Since the submission to the PDB, the docking predictions were
docking or assembly modeling performance may cru- shown to provide useful information, which often con-
cially depend on the accuracy of the models of individual firmed the final assignment or helped resolve ambiguous
subunits, the targets chosen for this experiment were ones. This occurred for both homodimer and homote-
proteins deemed to be readily modeled using templates tramer targets.
from the PDB. Interestingly, templates were used mainly Lastly, we find that the docking prediction perform-
to model the structures of individual subunits, to limit ance for the genuine homodimer targets was superior to
the sampling space of docking solution or to filter such that obtained for heterocomplexes in previous CAPRI
these solutions. Only a few groups carried out template- rounds, in line with the expectation that, owing to their
based docking for the majority of the targets, and two of higher binding affinity (and larger and more hydropho-
those ranked amongst the top performers, indicating that bic interfaces), homodimers are easier to predict than
this relatively recent modeling strategy has potential. heterodimers. Much poorer prediction performance was
As part of our assessment we established that the accu- however achieved for genuine tetrameric targets, where
racy of the models of the individual subunits was an the inaccuracy of the homology-built subunit models
important factor contributing to high accuracy predic- and the smaller pair-wise interfaces limited the predic-
tions of the corresponding complexes. At the same time tion performance. Accurately modeling of higher order
we observed that highly accurate models of the protein assemblies from sequence information is thus an area
components are not necessarily required for identifying where progress is needed.
their association modes with acceptable accuracy.
Furthermore, we provide evidence that protein dock- ACKNOWLEDGMENTS
ing procedures and in some cases also specialized
We are most grateful to the PDBe at the European
template-based methods generally outperform off-the-
Bioinformatics Institute in Hinxton, UK, for hosting the
shelf template-based prediction of complexes. These find-
CAPRI website. Our deepest thanks go to all the struc-
ings apply to templates identified on the basis of
tural biologists and to the following structural genomics
sequence information alone, as well as to templates
initiatives: Northeast Structural Genomics Consortium,
structurally more similar to the target. The added value
Joint Center for Structural Genomics, NatPro PSI:Biol-
of docking methods was particularly significant for the
ogy, New York Structural Genomics Research Center,
more difficult targets, where the structures of the identi-
Midwest Center for Structural Genomics, Structural
fied best templates differed more significantly from the
Genomics Consortium, for contributing the targets for
target structure
this joint CASP-CAPRI experiment. MFL acknowledges
Thus, the assessment results presented here confirm
support from the FRABio FR3688 Research Federation
that the prediction of homodimer assemblies by homol-
“Structural & Functional Biochemistry of Biomolecular
ogy modeling techniques and docking calculations is fea-
Assemblies.”
sible, especially for stable dimers that feature interface
areas of 1000–1500 A˚2, whose size is comparable or
REFERENCES
larger than the one associated with transient heterocom-
plexes. They also confirm that docking procedures can
1.Alberts B. The cell as a collection of protein machines: preparing
represent a competitive advantage over standard homol-
thenextgenerationofmolecularbiologists.Cell1998;92:291–294.
ogy modeling techniques, when those are applied with- 2.Berman HM, Battistuz T, Bhat TN, Bluhm WF, Bourne PE,
out further improvements to model the complex. Burkhardt K, Feng Z, GillilandGL, Iype L, Jain S,Fagan P, Marvin
On the other hand, difficulties arise when the subunit J, Padilla D, Ravichandran V, Schneider B, Thanki N, Weissig H,
WestbrookJD,ZardeckiC.TheProteinDataBank.ActaCrystallogr
interface in the target is similar in size to those associ-
SectionDBiolCrystallogr2002;58:899–907.
45
ated with crystal contacts. Such cases were associated
3.Smith MT, Rubinstein JL. Structural biology. Beyond blob-ology.
with a number of targets where the oligomeric state Science2014;345:617–619.
346 PROTEINS

PredictionofHomoandHeteroproteinComplexesbyProteinDockingandModeling
4.Kundrotas PJ, Zhu Z, Janin J, Vakser IA. Templates are available to Grudinin S, Derevyanko G, Mitchell JC, Wieting J, Kanamori E,
model nearly all complexes of structurally characterized proteins. Tsuchiya Y, Murakami Y, Sarmiento J, Standley DM, Shirota M,
ProcNatlAcadSciUSA2012;109:9438–9441. Kinoshita K, Nakamura H, Chavent M, Ritchie DW, Park H, Ko J,
5.Berman HM, Coimbatore Narayanan B, Di Costanzo L, Dutta S, Lee H, Seok C, Shen Y, Kozakov D, Vajda S, Kundrotas PJ, Vakser
Ghosh S, Hudson BP, Lawson CL, Peisach E, Prlic A, Rose PW, IA, Pierce BG, Hwang H, Vreven T, Weng Z, Buch I, Farkash E,
ShaoC,YangH, Young J, Zardecki C.Trendspottingin the Protein Wolfson HJ, Zacharias M, Qin S, Zhou HX, Huang SY, Zou X,
DataBank.FEBSLett2013;587:1036–1045. WojdylaJA,KleanthousC,WodakSJ.Blindpredictionofinterfacial
6.Marks DS, Hopf TA, Sander C. Protein structure prediction from waterpositionsinCAPRI.Proteins2014;82:620–632.
sequencevariation.NatBiotechnol2012;30:1072–1080. 18.Lensink MF, Mendez R, Wodak SJ. Docking and scoring protein
7.OvchinnikovS,Kamisetty H,Baker D. Robustand accurate predic- complexes:CAPRI3rdEdition.Proteins2007;69:704–718.
tion of residue-residue interactions across protein interfaces using 19.Lensink MF, Wodak SJ. Docking and scoring protein interactions:
evolutionaryinformation.eLife2014;3:e02030 CAPRI2009.Proteins2010;78:3073–3084.
8.Whitehead TA, Baker D, Fleishman SJ. Computational design of 20.Goodsell DS, Olson AJ. Structural symmetry and protein function.
novel protein binders and experimental affinity maturation. Meth- AnnuRevBiophysBiomolStruct2000;29:105–153.
odsEnzymol2013;523:1–19. 21.Kuhner S, van Noort V, Betts MJ, Leo-Macias A, Batisse C, Rode
9.Whitehead TA, Chevalier A, Song Y, Dreyfus C, Fleishman SJ, De M,YamadaT,MaierT,BaderS,Beltran-AlvarezP,Castano-DiezD,
Mattos C, Myers CA, Kamisetty H, Blair P, Wilson IA, Baker D. Chen WH, Devos D, Guell M, Norambuena T, Racke I, Rybin V,
Optimization of affinity, specificity and function of designed influ- SchmidtA,YusE,AebersoldR,HerrmannR,BottcherB,Frangakis
enzainhibitorsusingdeepsequencing.NatBiotechnol2012;30:543–
AS, Russell RB, Serrano L, Bork P, Gavin AC. Proteome organiza-
548.
tioninagenome-reducedbacterium.Science2009;326:1235–1240.
10.Ward AB, Sali A, Wilson IA. Biochemistry. Integrative structural
22.Hura GL, Menon AL, Hammel M, Rambo RP, Poole FL, 2nd
biology.Science2013;339:913–915.
Tsutakawa SE, Jenney FE, Jr., Classen S, Frankel KA, Hopkins RC,
11.Wodak SJ, Janin J. Structural basis of macromolecular recognition.
YangSJ,ScottJW,DillardBD,AdamsMWTainerJA.Robust,high-
AdvProteinChem2002;61:9–73.
throughputsolutionstructuralanalysesbysmallangleX-rayscatter-
12.Ritchie DW. Recent progress and future directions in protein-
ing(SAXS).NatMethods2009;6:606–612.
proteindocking.CurrProteinPeptSci2008;9:1–15.
23.Krissinel E, Henrick K. Inference of macromolecular assemblies
13.Vajda S, Kozakov D. Convergence combination of methods in
fromcrystallinestate.JMolBiol2007;372:774–797.
protein-proteindocking.CurrOpinStructBiol2009;19:164–170.
24.Webb B, Sali A. Protein structure modeling with MODELLER.
14.Fleishman SJ, Whitehead TA, Strauch EM, Corn JE, Qin S, Zhou
MethodsMolBiol2014;1137:1–15.
HX, Mitchell JC, Demerdash ON, Takeda-Shitaka M, Terashi G,
25.Arnold K, Bordoli L, Kopp J, Schwede T. The SWISS-MODEL
MoalIH,Li X, BatesPA, ZachariasM, Park H,KoJS, LeeH,Seok
workspace: a web-based environment for protein structure homol-
C, Bourquard T, Bernauer J, Poupon A, Aze J, Soner S, Ovali SK,
ogymodelling.Bioinformatics2006;22:195–201.
OzbekP,TalNB,HalilogluT,HwangH,VrevenT,PierceBG,Weng
26.Song Y, DiMaio F, Wang RY, Kim D, Miles C, Brunette T,
Z, Perez-Cano L, Pons C, Fernandez-Recio J, Jiang F, Yang F, Gong
Thompson J, Baker D. High-resolution comparative modeling with
X, Cao L, Xu X, Liu B, Wang P, Li C, Wang C, Robert CH,
RosettaCM.Structure2013;21:1735–1742.
GuharoyM,LiuS,HuangY,LiL,GuoD,ChenY,XiaoY,London
27.de Vries SJ, van Dijk M, Bonvin AM. The HADDOCKweb server
N, Itzhaki Z, Schueler-Furman O, Inbar Y, Potapov V, Cohen M,
fordata-drivenbiomoleculardocking.NatProtoc2010;5:883–897.
Schreiber G, Tsuchiya Y, Kanamori E, Standley DM, Nakamura H,
28.Macindoe G, Mavridis L, Venkatraman V, Devignes MD, Ritchie
Kinoshita K, Driggers CM, Hall RG, Morgan JL, Hsu VL, Zhan J,
DW. HexServer: an FFT-based protein docking server powered by
Yang Y, Zhou Y, Kastritis PL, Bonvin AM, Zhang W, Camacho CJ,
graphicsprocessors.NucleicAcidsRes2010;38:W445–449.
Kilambi KP, Sircar A, Gray JJ, Ohue M, Uchikoga N, Matsuzaki Y,
29.Pierce BG, Wiehe K, Hwang H, Kim BH, Vreven T, Weng Z.
Ishida T, Akiyama Y, Khashan R, Bush S, Fouches D, Tropsha A,
ZDOCK server: interactive docking prediction of protein-protein
Esquivel-Rodriguez J, Kihara D, Stranges PB, Jacak R, Kuhlman B,
complexes and symmetric multimers. Bioinformatics 2014;30:1771–
Huang SY, Zou X, Wodak SJ, Janin J, Baker D. Community-wide
1773.
assessment of protein-interface modeling suggests improvements to
30.Lyskov S, Gray JJ. The RosettaDock server for local protein-protein
designmethodology.JMolBiol2011;414:289–302.
15.Moretti R, Fleishman SJ, Agius R, Torchala M, Bates PA, Kastritis docking.NucleicAcidsRes2008;36:W233–238.
PL,Rodrigues JP, TrelletM, BonvinAM, CuiM, Rooman M, Gillis 31.Pierce B, Tong W, Weng Z. M-ZDOCK: a grid-based approach for
D,DehouckY,MoalI,Romero-DuranaM,Perez-CanoL,PallaraC, Cn symmetric multimer docking. Bioinformatics 2005;21:1472–
JimenezB,Fernandez-RecioJ,FloresS,PacellaM,PraneethKilambi 1478.
K, Gray JJ, Popov P, Grudinin S, Esquivel-Rodriguez J, Kihara D, 32.Szilagyi A, Zhang Y. Template-based structure modeling of protein-
Zhao N, Korkin D, Zhu X, Demerdash ON, Mitchell JC, Kanamori proteininteractions.CurrOpinStructBiol2014;24:10–23.
E, Tsuchiya Y, Nakamura H, Lee H, Park H, Seok C, Sarmiento J, 33.Kallberg M, Wang H, Wang S, Peng J, Wang Z, Lu H, Xu J. Tem-
Liang S, Teraguchi S, Standley DM, Shimoyama H, Terashi G, plate-based protein structure modeling using the RaptorX web
Takeda-Shitaka M, Iwadate M, Umeyama H, Beglov D, Hall DR, server.NatProtoc2012;7:1511–1522.
Kozakov D, Vajda S, Pierce BG, Hwang H, Vreven T, Weng Z, 34.Tuncbag N, Keskin O, Nussinov R, Gursoy A. Fast and accurate
Huang Y, Li H, Yang X, Ji X, Liu S, Xiao Y, Zacharias M, Qin S, modeling of protein-protein interactions by combining template-
Zhou HX, Huang SY, Zou X, Velankar S, Janin J, Wodak SJ, Baker interface-based docking with flexible refinement. Proteins 2012;80:
D.Community-wideevaluationofmethodsforpredictingtheeffect 1239–1249.
of mutations on protein-protein interactions. Proteins 2013;81: 35.Zemla A. LGA: A method for finding 3D similarities in protein
1980–1987. structures.NucleicAcidsRes2003;31:3370–3374.
16.LensinkMF, Wodak SJ. Docking, scoring, andaffinity predictionin 36.KryshtafovychA,MonastyrskyyB,FidelisK.CASPpredictioncenter
CAPRI.Proteins2013;81:2082–2095. infrastructure and evaluation measures in CASP10 and CASP
17.LensinkMF,MoalIH,BatesPA,KastritisPL,MelquiondAS,Karaca ROLL.Proteins2014;82Suppl2:7–13.
E, Schmitz C, van Dijk M, Bonvin AM, Eisenstein M, Jimenez- 37.Cozzetto D, Kryshtafovych A, Fidelis K, Moult J, Rost B,
Garcia B, Grosdidier S, Solernou A, Perez-Cano L, Pallara C, TramontanoA.Evaluationoftemplate-basedmodelsinCASP8with
Fernandez-Recio J, Xu J, Muthu P, Praneeth Kilambi K, Gray JJ, standardmeasures.Proteins2009;77Suppl9:18–28.
PROTEINS 347

M.F.Lensinketal.
38.ZemlaA,Venclovas,MoultJFidelisK.Processingandevaluationof 44.Mukherjee S, Zhang Y. MM-align: a quick algorithm for aligning
predictionsinCASP4.Proteins2001;Suppl5:13–21. multiple-chain protein complex structures using iterative dynamic
39.Remmert M, Biegert A, Hauser A, Soding J. HHblits: lightning-fast programming.NucleicAcidsRes2009;37:e83
iterativeproteinsequencesearchingbyHMM-HMMalignment.Nat 45.Bahadur RP, Chakrabarti P, Rodier F, Janin J. A dissection of spe-
Methods2012;9:173–175. cific and non-specific protein-protein interfaces. J Mol Biol 2004;
40.SodingJ. ProteinhomologydetectionbyHMM-HMMcomparison. 336:943–955.
Bioinformatics2005;21:951–960. 46.Dey S, Pal A, Chakrabarti P, Janin J. The subunit interfaces of
41.Hildebrand A, Remmert M, Biegert A, Soding J. Fast and accurate weakly associated homodimeric proteins. J Mol Biol 2010;398:146–
automatic structure prediction with HHpred. Proteins 2009;77 160.
Suppl9:128–132. 47.Mackinnon SS, Malevanets A, Wodak SJ. Intertwined associations
42.Sali A, Blundell TL. Comparative protein modelling by satisfaction in structures of homooligomeric proteins. Structure 2013;21:638–
ofspatialrestraints.JMolBiol1993;234:779–815. 649.
43.Zhang Y, Skolnick J. TM-align: a protein structure alignment algo- 48.Lensink MF, Wodak SJ. Blind predictions of protein interfaces by
rithmbasedontheTM-score.NucleicAcidsRes2005;33:2302–2309. dockingcalculationsinCAPRI.Proteins2010;78:3085–3095.
348 PROTEINS