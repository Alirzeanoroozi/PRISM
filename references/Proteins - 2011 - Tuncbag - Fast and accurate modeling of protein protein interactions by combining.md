proteins
STRUCTUREOFUNCTIONOBIOINFORMATICS
Fast and accurate modeling of
protein–protein interactions by
combining template-interface-based
docking with flexible refinement
Nurcan Tuncbag,1 Ozlem Keskin,1* Ruth Nussinov,2,3 and Attila Gursoy1*
1
CenterforComputationalBiologyandBioinformatics,CollegeofEngineering,KocUniversity,34450Sariyer,
Istanbul,Turkey
2
BasicScienceProgram,SAIC-Frederick,Inc.,CenterforCancerResearchNanobiologyProgram,
NCI-Frederick,Frederick,Maryland21702
3
DepartmentofHumanGeneticsandMolecularMedicine,SacklerInstituteofMolecularMedicine,
SacklerSchoolofMedicine,TelAvivUniversity,TelAviv69978,Israel
ABSTRACT INTRODUCTION
Thesimilaritybetweenfoldingandbindingledustoposit Identification of protein–protein interactions at the struc-
the concept that the number of protein–protein interface tural level is pivotal to the understanding of protein function.
motifs in nature is limited, and interacting protein pairs Experimentally, X-ray crystallography and NMR techniques
can use similar interface architectures repeatedly, even if are used to determine structures at atomic resolution. Despite
their global folds completely vary. Thus, known protein– the fast growth in the number of experimentally solved struc-
protein interface architectures can be used to model the tures of protein complexes, they are vastly fewer when com-
complexes between two target proteins on the proteome
pared to the available pair-wise interactions coming from
scale, even if their global structures differ. This powerful
high-throughput techniques. For this reason, computationally
conceptiscombinedwithaflexiblerefinementandglobal
fast and accurate algorithms to model protein complexes are
energy assessment tool. The accuracy of the method is
1,2
increasingly needed. The details of protein interactions pro-
highly dependent on the structural diversity of the inter-
vide crucial information, which allows fuller representation of
face architectures in the template dataset. Here, we vali-
date this knowledge-based combinatorial method on the the cellular network and the design of therapeutic molecules.
Docking Benchmark and show that it efficiently finds Computationally, structural modeling of protein interac-
high-quality models for benchmark complexes and their tions is performed predominantly in two ways: docking and
bindingregionsevenintheabsenceoftemplateinterfaces knowledge-based prediction (reviewed in Ref. 3). Docking is
having sequence similarity to the targets. Compared to the most widely used approach to predict the interactions of a
‘‘classical’’ docking, it is computationally faster; as the given pair of unbound protein structures. 4–9 Despite the
number of target proteins increases, the difference 10,11
improvements in blind docking tests, scoring functions
becomes more dramatic. Further, it is able to distinguish are not yet optimized and still constitute a major hurdle.12,13
bindersfromnonbinders.Thesefeaturesallowperforming
Docking is challenging on large scale because of two main rea-
large-scalenetworkmodeling.Theresultsonanindepend-
sons. First and most important reason is that docking will
ent target set (proteins in the p53 molecular interaction
practically always find ‘‘reasonable’’ solutions with apparent
map) show that current method can be used to predict
whetheragivenproteinpairinteracts.Overall,whilecon- favorable interactions and shape complementarity. For this
strainedbythediversityofthetemplateset,thisapproach reason, for reliable docking, it is desirable to not only have
efficiently produces high-quality models of protein– experimental information that the pair of proteins interact
protein complexes. We expect that with the growing
number of known interface architectures, this type of
Grantsponsor:TUBITAK;Grantnumbers:109T343and109E207;Grantsponsor:Turkish
knowledge-basedmethodswillbeincreasinglyusedbythe
AcademyofSciences;Grantsponsors:IntramuralResearchProgramoftheNIH,National
broadproteomicscommunity. CancerInstitute,CenterforCancerResearch,FederalFundsfromNationalCancerInsti-
tute,NationalInstitutesofHealth;Grantnumber:HHSN261200800001E.
Proteins2012;80:1239–1249. *Correspondenceto:OzlemKeskinorAttilaGursoy,CenterforComputationalBiologyand
VVC 2011WileyPeriodicals,Inc. Bioinformatics,CollegeofEngineering,KocUniversity,RumelifeneriYolu,34450Sariyer,
Istanbul,Turkey.E-mail:okeskin@ku.edu.troragursoy@ku.edu.tr.
Received24July2011;Revised29November2011;Accepted13December2011
Keywords: template-based docking; flexible refinement;
Published online 22 December 2011 in Wiley Online Library (wileyonlinelibrary.com).
3D modeling; protein interaction prediction. DOI:10.1002/prot.24022
VVC 2011WILEYPERIODICALS,INC. PROTEINS 1239

 10970134, 2012, 4, Downloaded from https://onlinelibrary.wiley.com/doi/10.1002/prot.24022 by Koc University, Wiley Online Library on [11/02/2026]. See the Terms and Conditions (https://onlinelibrary.wiley.com/terms-and-conditions) on Wiley Online Library for rules of use; OA articles are governed by the applicable Creative Commons License
N.Tuncbagetal.
but also some biochemical data related to the location of flexible refinement using a new efficient docking
the binding site. Second, large-scale docking is computa- method.35 The energy assessments make the prediction
tionally very demanding. Currently, a main advantage of more physical and provide a way to score the modeled
docking is the refinement of the modeled complexes and complexes. This approach is being released as a protocol
assessment of their energies. to be used freely community wide. 34
Knowledge-based modeling of protein interactions uses In this work, we present the validation of this tem-
information related to known structures of protein–pro- plate-based docking method integrated with flexible
tein complexes to target protein structures whose interac- refinement (PRISM—prediction of protein interactions
tions are unknown. Such a strategy is based on the con- by structural matching) on the Docking Benchmark. As
cept that binding and folding are similar events, and pro- long as structurally similar template interfaces are avail-
14–19
tein–protein interactions resemble protein cores. able in the dataset, it efficiently produces the near native
Similar to protein cores, the number of interface motifs or acceptable protein complex models independent from
is also limited in nature. 14,15,19–21 This type of model- the global folds. Further, unlike docking, it can also be
ing started with the pioneering work of Aloy and Russell used for prediction of pair-wise protein interactions,
where global structural homology of the templates was besides modeling binding pose and prediction of binding
22
considered in the modeling. Several works that use site. To show this feature, verification of the predicted
global sequence or structural similarities have also pair-wise interactions on large scale is performed on the
appeared (reviewed in Ref. 3). Later, the concept that proteins compiled and organized in the p53 molecular
36
binding sites of the proteins are structurally more con- interaction map (MIM). This template-based method
served than the remaining surface regions 23–25 has also is computationally very fast compared to classical dock-
been introduced. Further, it has been shown that unique ing on large scale. The running time analysis shows that
interface architectures can be used repeatedly by different its computational effectiveness is getting more obvious
protein pairs independent from their global folds 14,26 with an increasing number of target proteins.
| (reviewed        | in Ref. 27). | These | theoretical   | considerations |        |     |     |     |     |
| ---------------- | ------------ | ----- | ------------- | -------------- | ------ | --- | --- | --- | --- |
| and observations | motivated    |       | the idea that | interface      | archi- |     |     |     |     |
METHODS
| tectures         | can be utilized | instead | of global          | similarity | to        |     |     |     |     |
| ---------------- | --------------- | ------- | ------------------ | ---------- | --------- | --- | --- | --- | --- |
| infer structural | models          | of      | protein complexes. |            | The first |     |     |     |     |
Targetdataset
| approach, | based on | interface | templates, | was released | by        |     |     |     |     |
| --------- | -------- | --------- | ---------- | ------------ | --------- | --- | --- | --- | --- |
| et        | al. 28   |           |            |              | et al. 29 |     |     |     |     |
Aytuna (and its webserver by Ogmen ); Eighty-eight rigid-body test cases (from 165 protein
there, two complementary parts of the template-interface chains) in Docking Benchmark 3.0 37 are used for valida-
| architectures | served as | templates | in search | for structurally |     |                     |               |          |            |
| ------------- | --------- | --------- | --------- | ---------------- | --- | ------------------- | ------------- | -------- | ---------- |
|               |           |           |           |                  |     | tion of the method. | The benchmark | contains | 28 enzyme/ |
similar target protein surfaces. This approach handles inhibitor, 21 antibody/antigen, and 39 other type of com-
both spatially local and implicitly global similarities. plexes. All possible pairs of the 165 target protein chains
|     |     |     | appeared30,31 |     |     |     |     | 3   |     |
| --- | --- | --- | ------------- | --- | --- | --- | --- | --- | --- |
Later, additional similar approaches with are searched on the templates, and a 165 165 interac-
some differences, such as the size and type of the tem- tion matrix is constructed to see, if the method distin-
| plate dataset, | and the | structural | alignment | method. |     | A               |                  |     |     |
| -------------- | ------- | ---------- | --------- | ------- | --- | --------------- | ---------------- | --- | --- |
|                |         |            |           |         |     | guishes binders | from nonbinders. |     |     |
recent study showed that an interface architecture can To show the pathway-scale performance of the method,
also be used to perform multiple functions, 32 which sup- theproteinstructuresavailableintheMIM 36 isused.Sev-
ports the template-based modeling approaches. Further, eral proteins do not have complete structures, only frag-
these approaches have been found to be sufficiently accu- ments. For example, the full-length human DNA excision
| rate for | proteome scale | studies.33 | Another | advantage |     | as            |                   |          |               |
| -------- | -------------- | ---------- | ------- | --------- | --- | ------------- | ----------------- | -------- | ------------- |
|          |                |            |         |           |     | protein ERCC1 | has 297 residues; | however, | the available |
compared to docking is that template-based approaches structures are for residues 96–227 (PDB: 2a1i, chain A)
can predict whether two given proteins interact (or not). and 220–297 (1z00, chain A). Both fragments are consid-
Nonetheless, in all cases, the predicted interactions need ered in the target set. In this pathway, 77 proteins have
refinements, because all template-based methods consider structural information but when considering all protein
| only rigid-body | structural | alignment | for | docking | of the |     |     |     |     |
| --------------- | ---------- | --------- | --- | ------- | ------ | --- | --- | --- | --- |
fragments,thenumberofchainsincreasesto112.
targetproteins.Proteinsarenotrigidmolecules;theorien-
| tation of | the side chains | and | the movements | of  | the back- |     |     |     |     |
| --------- | --------------- | --- | ------------- | --- | --------- | --- | --- | --- | --- |
Templatedataset
| bone need | to be taken | into | account during | or  | following |     |     |     |     |
| --------- | ----------- | ---- | -------------- | --- | --------- | --- | --- | --- | --- |
therigid-bodyalignment. Three template sets are used throughout the work.
A multiscale strategy, which combines knowledge- Only the first and second are used for benchmarking: (i)
based structural alignments and docking, can lead to a An optimal template set for target proteins in the Dock-
34
powerful method to model the structural proteome. ing Benchmark: this set is extracted from the bound
The predicted complexes in which the interacting states of the proteins resulting in 88 interfaces, which
proteins have surfaces matching experimental interfaces contain discontinuous residue subsets of the target pro-
in a structurally nonredundant template dataset undergo tein chains. Using benchmark templates, we first check if
1240 PROTEINS

 10970134, 2012, 4, Downloaded from https://onlinelibrary.wiley.com/doi/10.1002/prot.24022 by Koc University, Wiley Online Library on [11/02/2026]. See the Terms and Conditions (https://onlinelibrary.wiley.com/terms-and-conditions) on Wiley Online Library for rules of use; OA articles are governed by the applicable Creative Commons License
Template-BasedDockingandFlexibleRefinement
the method can find all ‘‘true’’ solutions and no ‘‘false’’ appropriate for structural alignment of target surfaces to
ones with the optimal template set. (ii) A more diverse template partners. Geometry and residue type (hydro-
template dataset, which is utilized for the validation of phobic, hydrophilic, aromatic, or glycine) are considered
the method on the benchmark proteins: this second tem- in the structural alignment. Forty percent of the residues
plate set is constructed from a nonredundant interface of template chains should geometrically match the target
dataset composed of 49,512 interfaces structurally surfaces to pass to the next step. This threshold is 60%
clustered into 8205 different interface architectures. Elim- for template chains containing less than 50 residues. To
inating nonprotein complexes from these 8205 protein get rid of interfaces, which have small contact areas that
interfaces resulted in 7922 protein–protein interfaces. We could reflect crystal packing effects or spurious match-
aim to see how many interactions the method predicts ings, at least 15 residues should be matched in both
with this unbiased template set for our targets. Also, we cases. If there are computational hot spots in the tem-
selected heterodimeric protein interfaces (1036 interfaces) plate interface, at least one hot spot in each template
among this dataset to be used for the MIM analysis. (iii) partner should correctly match with the target surface.
A template set, composed of the interfaces already avail- Hot spot filtering incorporates evolutionary similarity
able in the MIM (59 interfaces): this set is used for the between the target surface and template interface in addi-
38
verification of other interactions in MIM. tion to structural similarity. The maximal RMSD
A˚.
Here, the interface definition is as follows: the contacts allowed in the structural alignment is 2 Thus, while
between two complementary chains are calculated from the alignment is rigid, it can handle movements within
A˚.
the distance between any two atoms each from one 2 Candidate target proteins passing the alignment
chain. If this distance is less than a threshold, which is threshold are transformed onto the corresponding tem-
defined as the sum of van der Waals radii of the corre- plate interface. If the two partners present (more than
A˚,
sponding atoms plus 0.5 they are considered as con- five) spatially colliding residues after tranformation, the
tacting residues. The neighbors of the contacting residues match is eliminated. Side-chain clashes are not consid-
| are searched | within the | same | chain and | the | threshold for | ered at this | stage. |     |     |     |     |
| ------------ | ---------- | ---- | --------- | --- | ------------- | ------------ | ------ | --- | --- | --- | --- |
Ca
| the distance    | between     | the | atoms is 6    | A˚. These | residues, |     |     |     |     |     |     |
| --------------- | ----------- | --- | ------------- | --------- | --------- | --- | --- | --- | --- | --- | --- |
| called ‘‘nearby | residues,’’ |     | are important | for       | a correct |     |     |     |     |     |     |
Flexiblerefinement
| structural      | alignment | of  | the interface |             | architec- |          |            |     |           |         |           |
| --------------- | --------- | --- | ------------- | ----------- | --------- | -------- | ---------- | --- | --------- | ------- | --------- |
| tures. 14,19,26 | Hot spots | in  | the template  | interfaces, | that      |          |            |     |           |         |           |
|                 |           |     |               |             |           | Flexible | refinement | of  | the rigid | docking | solutions |
is, residues contributing more to the binding energy, are involves resolving steric clashes, especially of side chains,
38
identified using the Hotpoint web server. followed by ranking putative complexes by the global
35
|     |     |     |     |     |     | energy.         | FiberDock | is used    | for | flexible | refinement,  |
| --- | --- | --- | --- | --- | --- | --------------- | --------- | ---------- | --- | -------- | ------------ |
|     |     |     |     |     |     | which considers | both      | side-chain | and | backbone | flexibility. |
Thepredictionalgorithm
|     |     |     |     |     |     | The side-chain | orientations |     | are | optimized | using a |
| --- | --- | --- | --- | --- | --- | -------------- | ------------ | --- | --- | --------- | ------- |
Asnotedabove,thisalgorithmcombinestemplate-based rotamer library, and the combination of rotamers having
predictionwithflexiblerefinement.Proteins interactusing lowest total energy is selected. Side-chain optimization is
| their surfaces. | Different | from | other approaches, |     | instead of |            |         |              |          |     |                  |
| --------------- | --------- | ---- | ----------------- | --- | ---------- | ---------- | ------- | ------------ | -------- | --- | ---------------- |
|                 |           |      |                   |     |            | restricted | to only | the clashing | residues |     | in the predicted |
considering the overall structure, we extract the surface interface. Up to 20% clashes between side chains are
residuesofthetarget proteins.Thispreventsaninaccurate allowed. Backbone flexibility is modeled using the first
matching of the templates to the core of the proteins, 50 normal modes of the corresponding protein. The
especially in proteins with large sizes. Here, we define the quality of the predicted models is assessed using the
| surface | as a shell around | the | protein. The | surface | residues |            |               |        |       |     |                  |
| ------- | ----------------- | --- | ------------ | ------- | -------- | ---------- | ------------- | ------ | ----- | --- | ---------------- |
|         |                   |     |              |         |          | calculated | global energy | value, | which | is  | the single score |
are found by calculating the relative accessibilities of the to rank the models. The lowest global energy value
residues using NACCESS.39 If the relative accessibility, implies the highest ranking prediction. Other thresholds
that is, the ratio of the accessible surface to the maximum for residue match ratio, hot spots matching, and geomet-
accessibilityofthatresidueinanextendedpeptideconfor- rical clashes in the previous steps are used only for filter-
| mation, | is more than | 15%,   | these residues       | are | defined    | as               |               |     |     |     |     |
| ------- | ------------ | ------ | -------------------- | --- | ---------- | ---------------- | ------------- | --- | --- | --- | --- |
|         |              |        |                      |     |            | ing the possible | interactions. |     |     |     |     |
| surface | residues. As | in the | template interfaces, |     | structural |                  |               |     |     |     |     |
scaffoldsonthesurfacesareveryimportantforanaccurate
matching;thus,nearbyresiduesofsurfaceresiduesarealso RESULTS AND DISCUSSIONS
calculated,asdescribedearlier.
PRISMeffectivelyfindsnearnativemodels
anddistinguishesthenonbindersfrom
| Rigid-bodyalignment |     |     |     |     |     | bindersonanoptimaltemplateset |     |     |     |     |     |
| ------------------- | --- | --- | --- | --- | --- | ----------------------------- | --- | --- | --- | --- | --- |
Structural aligment of surface regions requires a At the validation stage, our first aim is to examine
method that can compare discontinuous fragments in a how the method performs on an optimal template set.
sequence-order-independent fashion. 40,41 MultiProt 41 is As expected, this is the best scenario that the method can
|     |     |     |     |     |     |     |     |     |     | PROTEINS | 1241 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | -------- | ---- |

N.Tuncbagetal.
complexes. The network representation of these extra
interactions is illustrated in Figure 1. Most of these arise
from antibodies or antigens. Prediction of antibody/anti-
gen complexes is challenging; many efforts in library con-
struction aimed to understand the binding specificity of
42
antibodies to antibody.
As an example for the extra interactions, the modeled
complexes of bovine trypsin are illustrated. In addition
to the soybean trypsin inhibitor (1ba7), our algorithm
predicts that bovine trypsin (1qqu) can interact with
Bowman–Birk inhibitor (1k9b:A), pancreatic secretory
trypsin inhibitor (1hpt), bovine pancreatic trypsin inhibi-
tor (9pti), CMTI-1 squash (1lu0:B) inhibitor, and TDPI
from tick (2uux), which are all trypsin inhibitors.
Although the overall structures and sequences of the
partner proteins are dissimilar, they can bind to the
bovine trypsin on the same surface, and the energy-based
rankings of these interactions are high. In Figure 2, the
partners of the bovine trypsin are superimposed onto
Figure 1 each other to illustrate clearly the structural similarity in
Illustrationoftheextrainteractionspredictedbyourmethodonthe the interface region.
benchmarktemplates(88interfaces).The165proteinchainsare
categorizedintofiveclassesinlinewithDockingBenchmark
classification.Eachnoderepresentsoneclass,andeachedgerepresents
PRISMfindshigh-qualitymodelsforthe
thenumberofpredictedinteractionsbetweentwoclasses.[Colorfigure
proteinsintheDockingBenchmark
canbeviewedintheonlineissue,whichisavailableat
independentfromtheglobalfoldsofthe
wileyonlinelibrary.com.]
templateinterfaces
face, because all the templates are coming from the native We next use an unbiased and more diverse template
complexes of the target proteins. Although running the dataset containing 7922 protein interfaces, which are
method on templates coming from ‘‘self-hits’’ seems
redundant, it provides important information about the
performance, such as, learning how the structural align-
ment is performed and selecting the optimum parameters
for alignment, because both target surfaces and template
partners are discontinuous sets of residues. This analysis
also provides information on how the method distin-
guishes binders from nonbinders. Here, the aim is to ver-
ify the matching parameters and to show that given a
template set containing similar interfaces, the method
finds near native modes of the protein complexes with
relatively less false positives. Each of the target protein
surfaces is aligned with the partner chains of those 88
interfaces. The method is applied to all possible pairs
(165 3 165), and a matrix of interacting pairs is gener-
ated. At the matching phase, correct binding regions are
found for all 88 protein complexes, except one case,
which is an antibody/antigen complex (Fab N10/Staphy-
lococcal nuclease complex). These correct protein com-
plex models are highly ranked by FiberDock. Besides the
Figure 2
87 complex models, 243 protein complexes are also mod-
Bovinetrypsin(coloredwhite)caninteractwithseveraltrypsin
eled at the end of this run. If all 165 nodes would inter- inhibitorsusingthesameregionandthreeofthesepartnersare
act with each other, there would be 13,530 edges in the superimposedtoshowthestructuralsimilarityintheirbindingsites
network. Our algorithm gives only 243 extra interactions; only.AlthoughtheoverallstructuresofBowman–Birkinhibitor
(1k9b:A,yellow),bovinepancreatictrypsininhibitor(9pti,pink),and
41 of them are modeled as antibody/antigen, 55 as
TDPIfromtick(2uux,cyan)aredissimilar,thebindingregionto
enzyme and inhibitor/substrate complexes, 74 as one-side
bovinetrypsinisstructurallyveryconserved.[Colorfigurecanbe
antibody, and the remaining are between other types of viewedintheonlineissue,whichisavailableatwileyonlinelibrary.com.]
1242 PROTEINS
10970134,
2012,
4,
Downloaded
from
https://onlinelibrary.wiley.com/doi/10.1002/prot.24022
by
Koc
University,
Wiley
Online
Library
on
[11/02/2026].
See
the
Terms
and
Conditions
(https://onlinelibrary.wiley.com/terms-and-conditions)
on
Wiley
Online
Library
for
rules
of
use;
OA
articles
are
governed
by
the
applicable
Creative
Commons
License

Template-BasedDockingandFlexibleRefinement
Figure 3
DistributionofI-Scoreversusf andiRMSDvaluesbasedon(a)defaultparametersand(b)hotspotthresholdrelaxedparameters.HigherI-Score
nat
implieshigherfractionofnativecontacts.TheiRMSDvaluesofhighI-Scoresfluctuatebetween0and3.[Colorfigurecanbeviewedintheonline
issue,whichisavailableatwileyonlinelibrary.com.]
structurally nonredundant. We eliminated template inter- while in the relaxed case, it produces 136 binding modes
faces, if at least one side of these interfaces is similar in on average for each target pair. Therefore, the distribu-
sequence to one of the target proteins. tion in Figure 3(a) is less populated when compared to
The comparison with the native complex and the scor- Figure 3(b). This distribution shows that f and I-Score
nat
ing of the quality of the binding mode are performed by are linearly correlated with each other. Further, predicted
the classical docking scoring parameters, such as lRMSD complexes having high I-Score have low iRMSD values.
(the RMSD between the native and modeled ligands after As a result, measuring the quality of the predicted mod-
superimposition of the receptors), iRMSD (the RMSD els using I-Score performs well.
between the native and modeled interfaces), and f The template dataset is not biased to the native com-
nat
(fraction of the native contacts). I-Score, a new metric to plexes in the Docking Benchmark, but it includes several
43
score the quality of docking predictions, is also used to self-hits in it. We remove these from the template data-
compare the predicted model with the native complex. If set. In the default case, for 25 out of 88 targets pairs, at
the I-Score is greater than 0.17, the predicted binding least one near native conformation is found (I-Score
mode is a near native model, and if it is between 0.12 (cid:1)0.17); also two more complexes are predicted as ac-
and 0.17, it is an acceptable model. If I-Score is less than ceptable models (0.12 < I-Score < 0.17). In the relaxed
0.12, the predicted model is incorrect. We checked the version, 41 target pairs are modeled as near native
distribution of the I-Score versus f and I-Score versus (I-Score (cid:1) 0.17) after removal of self-hits. In addition,
nat
iRMSD values for the modeled complexes by PRISM at least one acceptable model is produced for 25 target
(Fig. 3). Here, we performed two distinct PRISM runs. pairs among the remaining ones (0.12 < I-Score <0.17).
The former is performed using default thresholds; the To see how our method performs on target proteins, if
latter is performed by relaxing the hot spot matching we eliminate homologous template chains, we tried
threshold. In the default case, PRISM produces seven different sequence similarity thresholds between target
binding modes on average for each target protein pair, proteins and template-interface partners. In Figure 4(a),
PROTEINS 1243
10970134,
2012,
4,
Downloaded
from
https://onlinelibrary.wiley.com/doi/10.1002/prot.24022
by
Koc
University,
Wiley
Online
Library
on
[11/02/2026].
See
the
Terms
and
Conditions
(https://onlinelibrary.wiley.com/terms-and-conditions)
on
Wiley
Online
Library
for
rules
of use;
OA
articles
are
governed
by
the
applicable
Creative
Commons
License

N.Tuncbagetal.
Figure 4
(a)ThechangeinthetotalnumberofnearnativeandacceptablemodelspredictedbyPRISMwiththedefaultparametersandthehotspot
thresholdrelaxedparametersversusthesequencesimilarityeliminationthresholdsbetweenthetargetproteinsandtemplate-interfacepartners.At
100%sequencesimilarity,templateinterfacesfromnativecomplexesintheDockingBenchmarkareexcluded.(b)ThedistributionoftheI-Score
versuscalculatedglobalenergies(DG )andI-Scoreversuscombinedmatchingscore(totalof0.63f 10.43f foreachmatchingside
calc hotspot match
ofthetemplateinterface)valuesforpredictedcomplexes.[Colorfigurecanbeviewedintheonlineissue,whichisavailableat
wileyonlinelibrary.com.]
the total number of near native and acceptable models The subtilisin (2gkr)/ovomucoid (1scn) complex is
predicted by PRISM are illustrated versus the sequence modeled using the template interface between subtilisin/
similarity threshold. In Figure 4(b), the distributions of chymotrypsin inhibitor 2 (2sniEI). The sequence similar-
the I-Score versus calculated global energies (DG ) and ity between ovomucoid and chymotrypsin inhibitor is
calc
I-Score versus combined matching score (total of 0.6 3 low (only 7%). Although their global folds are dissimilar,
f 1 0.4 3 f for each matching side of the tem- the structural similarity between the binding regions is
hotspot match
plate interface) values for the modeled complexes by very high. The global energy is calculated as 248.57 kcal/
PRISM are shown, where f is the fraction of mol for this interaction. To show the similarity between
hotspot
identically matched hot spot residues, and f is the the binding regions and dissimilarity in the global folds,
match
fraction of the matched residues. Also, predicted protein chymotrypsin inhibitor 2 is superimposed on ovomucoid,
complexes having negative energy values are mostly near and the predicted model is illustrated in Figure 5(a). The
native complexes. In addition, predicted protein com- iRMSD for this binding mode is 0.65 A˚, and the I-Score
plexes having a combined matching score greater than is 0.7585. Our method successfully identifies the binding
one are mostly near native. The matching scores and region on ovomucoid and correctly models the subtilisin/
calculated global energies are correlated with each other. ovomucoid complex.
Besides providing a ranking metric as global energy, The interaction between bovine chymotrypsinogen
flexible refinement solves side-chain clashes in the (2cga) and pancreatic secretory trypsin inhibitor (1hpt)
interface region, minimizes the overall protein complex; is found using the interface region in human leukocyte
in this way, produces physically meaningful models. elastase/the turkey ovomucoid inhibitor complex
Later, we will show several examples among the correctly (1ppfEI). The sequence similarity between elastase and
modeled protein complexes. chymotrypsinogen is 32% and between trypsin inhibitor
1244 PROTEINS
10970134,
2012,
4,
Downloaded
from
https://onlinelibrary.wiley.com/doi/10.1002/prot.24022
by
Koc
University,
Wiley
Online
Library
on
[11/02/2026].
See
the
Terms
and
Conditions
(https://onlinelibrary.wiley.com/terms-and-conditions)
on
Wiley
Online
Library
for
rules
of
use;
OA
articles
are
governed
by
the
applicable
Creative
Commons
License

Template-BasedDockingandFlexibleRefinement
Figure 5
(a)Thesubtilisin(2gkr,white)/ovomucoid(1scn,cyan)complexismodeledonthetemplatesubtilisin/chymotrypsininhibitor2(2sniEI).
Chymotrypsininhibitor2(pink)issuperimposedonovomucoidtoshowthestructuralsimilarityintheinterfaceregionbetweentargetand
templatechains.(b)Theinteractionbetweenbovinechymotrypsinogen(2cga,pink)andpancreaticsecretorytrypsininhibitor(1hpt,white)is
modeledontheinterfaceregionofhumanleukocyteelastase/theturkeyovomucoidinhibitorcomplex.Templateinterface(1ppfEI)iscoloredcyan
andgreentoshowthestructuralmatchingbetweentargetsurfacesandtemplatepartners.
and ovomucoid inhibitor is 28%. As illustrated in Figure matripase (1eax) and trypsin (2tpiZ) is 38% and between
5(b), the template-interface matches well with target matripase (1eax) and chymotrypsin (1t8oA) is 33%.
surfaces and the calculated global energy for this interac- Beta-trypsin and tryptase inhibitor complex is also
tion is 251.10 kcal/mol. The iRMSD value for this modeled on the template from thrombin–trypsin inhibi-
binding mode is 1.51 A˚, and the I-Score is 0.4918. This tor complex (1bthKQ). The sequence similarity between
complex is also predicted on the template interface in thrombin and beta-trypsin is 36% and between tryptase
rhodniin in complex with thrombin (1tbrHR). The inhibitor and trypsin inhibitor is 25%. This model is
sequence similarity between thrombin (1tbrH) and chy- identified as near native according to the scoring meas-
motrypsinogen (2cga) is 35% and between rhodniin and ures (lRMSD 5 2.66 A˚, iRMSD 5 0.62 A˚, I-Score 5
trypsin inhibitor is 25%. I-Score is 0.304 and shows that 0.62, and f 5 0.75). The calculated global energy for
nat
this predicted model is near native. The calculated global this model is 237.67 kcal/mol.
energy is 232.96 kcal/mol. The lRMSD is 2.73 A˚, and Although there is no sequence similarity, a high-quality
the iRMSD is 1.82 A˚ for this model. modeloftheribonuclease(1rghA)/barstar(1a19)complex
The template interface within idiotope–anti-idiotope is obtained on the template interface of catalytic antibody
complex (1cicBC) produced a near native model for the 4B2complex(1f3dHJ),wheretheI-Scoreis0.282,andthe
human CD40 ligand (1aly) and the immunoglobulin Fab iRMSD is 1.41 A˚. In Figure 6, the predicted complex is
fragment (1i9rH) complex. The sequence similarity superimposedontothenativecomplex.Thefractionofthe
between 1cicB and 1i9rH is 67% and between 1cicC and native contacts in the modeled complex is found to
1aly is 8%. The lRMSD measure is 7.6 A˚, and iRMSD is be 0.50. The calculated global energy for this model is
2.13 A˚ for this model. The I-Score is 0.329 and 58% of 218.25kcal/mol.Asshowninthefigure,whilethebinding
the native contacts are correctly predicted (f 5 0.58). region is highly overlapping with the native complex, the
nat
The calculated global energy is 13.18 kcal/mol. predictedorientationisshifted10.27A˚,whenwecompare
The matripase (1eax)/trypsin inhibitor (9pti) complex overallstructuresoftheligandproteins.
is modeled on the template interface extracted from Overall, the validation results show that even if we
thrombin/trypsin inhibitor (1bthKQ). The metrics for eliminateinterfaceshaving morethan 50%sequencesimi-
this model are as follows: I-Score 5 0.745, iRMSD 5 larityfromthetemplateset,themethodstillprovidesboth
0.58 A˚, and f 5 0.79. The sequence similarity between near native and acceptable models (17 out of 88 with
nat
thrombin and matripase is 35%. This complex is also defaultthresholdsand61outof88withhotspotthreshold
correctly modeled on the template interfaces within tryp- relaxed case) using this method independent from the
sinogen and trypsin inhibitor (2tpiZI, I-Score 5 0.71, sequence similarity and global folds of the corresponding
iRMSD 5 0.63 A˚, and f 5 0.79) and chymotrypsin template interfaces. The flexible refinement of these mod-
nat
and trypsin inhibitor (1t8oAB, I-Score 5 0.64, iRMSD 5 eled protein complexes makes these models physically
0.85 A˚, and f 5 0.72). The sequence similarity between more meaningful. We considered only rigid-body cases in
nat
PROTEINS 1245
10970134,
2012,
4,
Downloaded
from
https://onlinelibrary.wiley.com/doi/10.1002/prot.24022
by
Koc
University,
Wiley
Online
Library
on
[11/02/2026].
See
the
Terms
and
Conditions
(https://onlinelibrary.wiley.com/terms-and-conditions)
on
Wiley
Online
Library
for
rules
of
use;
OA
articles
are
governed
by
the
applicable
Creative
Commons
License

 10970134, 2012, 4, Downloaded from https://onlinelibrary.wiley.com/doi/10.1002/prot.24022 by Koc University, Wiley Online Library on [11/02/2026]. See the Terms and Conditions (https://onlinelibrary.wiley.com/terms-and-conditions) on Wiley Online Library for rules of use; OA articles are governed by the applicable Creative Commons License
N.Tuncbagetal.
|     |     |     |     |     |     |     | their corresponding |                      | protein    |          | structures, |            | their alignment | is          |
| --- | --- | --- | --- | --- | --- | --- | ------------------- | -------------------- | ---------- | -------- | ----------- | ---------- | --------------- | ----------- |
|     |     |     |     |     |     |     | faster than         | the                  | alignment  |          | of two      | complete   | proteins.       | On          |
|     |     |     |     |     |     |     | average,            | rigid-body           | structural |          | alignment   |            | of a            | target sur- |
|     |     |     |     |     |     |     | face and            | a template-interface |            |          | chain       | takes      | 1 s. To         | compare     |
|     |     |     |     |     |     |     | the first           | part                 | of this    | approach | with        | rigid-body |                 | docking     |
44
|     |     |     |     |     |     |     | algorithms, | we      | selected | PatchDock, |            |             | which       | performs   |
| --- | --- | --- | --- | --- | --- | --- | ----------- | ------- | -------- | ---------- | ---------- | ----------- | ----------- | ---------- |
|     |     |     |     |     |     |     | docking     | in less | than     | 10 min     | on average |             | on a single | proc-      |
|     |     |     |     |     |     |     | essor. In   | the     | running  | time       | analysis,  | we          | assume      | that the   |
|     |     |     |     |     |     |     | template    | set is  | composed |            | of 7922    | interfaces. |             | Also, each |
targettotemplatepartneralignmentisassumedtobeper-
|     |     |     |     |     |     |     | formed | in 1 s | as measured |     | earlier, | and | each docking | run |
| --- | --- | --- | --- | --- | --- | --- | ------ | ------ | ----------- | --- | -------- | --- | ------------ | --- |
foraproteinpairisassumedtobeperformedin10min.
|     |     |     |     |     |     |     | Docking     | must   | be          | performed | for       | all     | pairs  | of proteins |
| --- | --- | --- | --- | --- | --- | --- | ----------- | ------ | ----------- | --------- | --------- | ------- | ------ | ----------- |
|     |     |     |     |     |     |     |             |        | N           |           |           | N       | 3      | (N 2        |
|     |     |     |     |     |     |     | separately; | so     | for         | target    | proteins, |         |        | 1)/2        |
|     |     |     |     |     |     |     | docking     | runs   | are needed. | In        | total,    | docking | for    | N targets   |
|     |     |     |     |     |     |     |             |        | 3 N 3       | (N 2      |           |         | O(N2). |             |
|     |     |     |     |     |     |     | takes 10    | min    |             |           | 1)/2,     | which   | is     | On the      |
|     |     |     |     |     |     |     | other hand, | during | structural  |           | matching, |         | one    | target sur- |
Figure 6 face is compared to one side of a template interface only
|     |     |     |     |     |     |     |     |     |     |     |     |     | 3   | 3   |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
Thebarnase–barstarcomplexpredictedfromthetemplateinterface once. In total, structural matching takes 1s 2 1036
betweencatalyticantibody4B2complex.Thecyancoloredstructureis
|     |     |     |     |     |     |     | 3 N for | 1036 | templates | and | 1s  | 3 2 | 3 7922 | 3 N for |
| --- | --- | --- | --- | --- | --- | --- | ------- | ---- | --------- | --- | --- | --- | ------ | ------- |
barnase.Theyellow-coloredstructureispredictedorientationofthe
barstarandtransparentred-coloredoneisthenativeorientationofthe 7922 templates; that is, it increases linearly with increas-
O(N).
barstar.ThelRSMDvalueshowsthatthepredictedpartneris10A˚ ing number of targets, The rigid-body alignment
shiftedfromthenativestructure.Fiftypercentofthecontactsbetween part already defeats the docking in computational time.
barnase/barstararecorrectlyfoundandtheRMSDbetweentheinterface
ofpredictedandnativebarstar/barnasecomplexis1.41A˚.[Colorfigure Also, the total running time of the flexible refinement
|     |     |     |     |     |     |     | part is | dependent | on  | the | number | of  | output | solutions, |
| --- | --- | --- | --- | --- | --- | --- | ------- | --------- | --- | --- | ------ | --- | ------ | ---------- |
canbeviewedintheonlineissue,whichisavailableat
wileyonlinelibrary.com.] and the increase is again linear. Thus, even if we add the
|             |           |                |        |             |               |     | time spent  | for     | flexible   | refinement, |        | the          | difference | is still |
| ----------- | --------- | -------------- | ------ | ----------- | ------------- | --- | ----------- | ------- | ---------- | ----------- | ------ | ------------ | ---------- | -------- |
|             |           |                |        |             |               |     | very large. | Figure  | 7(a)       | illustrates |        | a comparison |            | of run-  |
| Docking     | Benchmark | 3.0            | during | validation, | because       | the |             |         |            |             |        |              |            |          |
|             |           |                |        |             |               |     | ning times  | as      | a function | of          | target | dataset      | size.      | As shown |
| first stage | of this   | template-based |        | method      | is rigid-body |     |             |         |            |             |        |              |            |          |
|             |           |                |        |             |               |     | in this     | figure, | for a      | small       | set of | target       | proteins,  | both     |
alignment.Therefore,ifatargetmonomerchangesitscon-
|           |               |                     |         |           |         |            | methods         | perform  | the         | docking      | in          | a rather  |      | equal time  |
| --------- | ------------- | ------------------- | ------- | --------- | ------- | ---------- | --------------- | -------- | ----------- | ------------ | ----------- | --------- | ---- | ----------- |
| formation | substantially | in                  | the     | binding   | region, | when       | it              |          |             |              |             |           |      |             |
|           |               |                     |         |           |         |            | frame. However, |          | at large    | scale,       | it changes, |           | and  | the knowl-  |
| binds to  | its partner   | proteins            | (called | difficult | cases), | it         | is              |          |             |              |             |           |      |             |
|           |               |                     |         |           |         |            | edge-based      | method   |             | dramatically |             | decreases | the  | solution    |
| hard to   | handle        | this conformational |         | change    | in      | the rigid- |                 |          |             |              |             |           |      |             |
|           |               |                     |         |           |         |            | space as        | a result | the running |              | times.      | As the    | size | of the tar- |
bodyalignmentstage.ThePRISMresultsondifficultcases
|                        |           |          |               |                |               |          | get dataset  | enlarges, | the        | difference |         | between    | running | times     |
| ---------------------- | --------- | -------- | ------------- | -------------- | ------------- | -------- | ------------ | --------- | ---------- | ---------- | ------- | ---------- | ------- | --------- |
| are consistent         | with      | this     | statement,    | with           | only          | two near |              |           |            |            |         |            |         |           |
|                        |           |          |               |                |               |          | gets larger. | Here,     | we         | show the   | time    | comparison |         | for up to |
| nativeandtwoacceptable |           |          | complexesthat |                | aremodeledout |          |              |           |            |            |         |            |         |           |
|                        |           |          |               |                |               |          | 165 target   | proteins  |            | in the     | Docking | Benchmark. |         | For a     |
| of the 17              | difficult | cases in | Docking       | Benchmark      |               | 3.0. As  | a            |           |            |            |         |            |         |           |
|                        |           |          |               |                |               |          | larger set,  | the       | difference | gets       | more    | dramatic.  |         | Also, the |
| future direction,      |           | we are   | working       | on integrating |               | flexible |              |           |            |            |         |            |         |           |
aligment into the first stage instead of rigid-body align- running time of our method is a function of the tem-
|     |     |     |     |     |     |     | plate dataset | size | in addition |     | to that | of  | the target | dataset. |
| --- | --- | --- | --- | --- | --- | --- | ------------- | ---- | ----------- | --- | ------- | --- | ---------- | -------- |
ment.Anotheroptioniscollectingalltheconformationsof
|     |     |     |     |     |     |     | In Figure | 7(b), | the | running | time | of  | the docking | for |
| --- | --- | --- | --- | --- | --- | --- | --------- | ----- | --- | ------- | ---- | --- | ----------- | --- |
thetemplatesavailableinthePDBandusingthemtohan-
|                    |     |         |         |     |         |     | 165 targets | is    | compared | to     | our      | template-based |     | method.    |
| ------------------ | --- | ------- | ------- | --- | ------- | --- | ----------- | ----- | -------- | ------ | -------- | -------------- | --- | ---------- |
| dle conformational |     | changes | between | the | unbound | and |             |       |          |        |          |                |     |            |
|                    |     |         |         |     |         |     | The figure  | shows | that     | if the | template | dataset        |     | size would |
boundstatesofthetargetproteins.
|     |     |     |     |     |     |     | be composed  |     | of 25,000 | interfaces, |       | the   | running | times of    |
| --- | --- | --- | --- | --- | --- | --- | ------------ | --- | --------- | ----------- | ----- | ----- | ------- | ----------- |
|     |     |     |     |     |     |     | both methods |     | would     | be the      | same, | which | is      | larger than |
Template-basedmodelingofprotein the current template sets. 14,19,21,26–28
complexesiscomputationallyfasterthan
|     |     |     |     |     |     |     | Our | method | is  | restricted | by  | the | diversity | of the |
| --- | --- | --- | --- | --- | --- | --- | --- | ------ | --- | ---------- | --- | --- | --------- | ------ |
‘‘classical’’docking
|     |     |     |     |     |     |     | interface | architectures. |     | Although |     | the number |     | of distinct |
| --- | --- | --- | --- | --- | --- | --- | --------- | -------------- | --- | -------- | --- | ---------- | --- | ----------- |
When modeling the interactions between target pro- interface architectures increases exponentially, not all
teins, most of the time is spent on rigid-body structural interface types have been discovered till date. As in all
alignment and flexible refinement of the filtered com- template-based methods (such as homology modeling,
plexes. Here, we show the running time differences motif finding-based predictions, and threading), the per-
between rigid-body alignment and rigid-body docking as formance of the method mainly depends on the quality
a function of the number of target proteins in the target of the template dataset. The contribution of the parame-
(N).
dataset Because the surface of the target protein and ters used during prediction loses its importance, when an
the interface of the template partners are only subsets of optimal template set is used. Overall, the running time
1246 PROTEINS

 10970134, 2012, 4, Downloaded from https://onlinelibrary.wiley.com/doi/10.1002/prot.24022 by Koc University, Wiley Online Library on [11/02/2026]. See the Terms and Conditions (https://onlinelibrary.wiley.com/terms-and-conditions) on Wiley Online Library for rules of use; OA articles are governed by the applicable Creative Commons License
Template-BasedDockingandFlexibleRefinement
|     |     | analysis     | shows           | the clear  | advantage  |             | of           | this method | in         |
| --- | --- | ------------ | --------------- | ---------- | ---------- | ----------- | ------------ | ----------- | ---------- |
|     |     | docking      | on a large      | scale.     | With       | the         | continuous   |             | growth of  |
|     |     | the PDB,     | knowledge-based |            | methods    |             | will         | get more    | attrac-    |
|     |     | tive for     | the community,  |            | especially |             | because      | of          | their fast |
|     |     | running      | times;          | further,   | the        | reliability |              | of using    | known      |
|     |     | motifs from  | nature          | will       | permit     | their       | application  |             | to pro-    |
|     |     | teins, which | a               | priori are | not        | known       | to interact. |             |            |
PRISMisapossibletoolforprediction
ofbindingpartnersofproteins:
apathway-scalecasestudyonp53MIM
|     |     | As above        | mentioned,     |                | we            | proposed | the     | predominance |            |
| --- | --- | --------------- | -------------- | -------------- | ------------- | -------- | ------- | ------------ | ---------- |
|     |     | of this method  |                | over classical |               | docking  | making  | it           | a possible |
|     |     | tool for        | the prediction |                | of whether    |          | any two | proteins     | inter-     |
|     |     | act. Therefore, |                | following      | validation,   |          | the     | next step    | is veri-   |
|     |     | fication        | of the         | pair-wise      | interactions. |          | For     | this         | purpose,   |
|     |     | we apply        | our multiscale |                | combinatorial |          | docking |              | algorithm  |
MIM.36
|          |     | to the              | proteins    | available      | in         | the         | human     |                | Addi-      |
| -------- | --- | ------------------- | ----------- | -------------- | ---------- | ----------- | --------- | -------------- | ---------- |
|          |     | tional interactions |             | from           | other      | databases   |           | such           | as DIP, 45 |
|          |     | 46                  |             | 47             |            | 48          |           |                |            |
|          |     | MINT,               | BIND,       | and            | IntAct     | enrich      | this      | map.           | Overall,   |
|          |     | 328 interactions    |             | between        | 104        | proteins    |           | are found      | from       |
|          |     | MIM and             | interaction |                | databases. | Among       |           | them,          | only for   |
|          |     | 25 interactions     |             | the structures |            | of the      | complexes |                | are avail- |
|          |     | able in             | the PDB.    | At this        | point,     | our         | method    | intervenes     | to         |
|          |     | complete            | the         | lacking        | network    | components. |           |                | Using the  |
|          |     | template            | interfaces  | in             | MIM        | with        | a default |                | matching   |
|          |     | threshold,          | 108         | interactions   |            | between     |           | 49 proteins    | are        |
|          |     | obtained,           | and         | 30 of          | these 108  | are         | known     | experimentally |            |
| Figure 7 |     | (without            | complex     | structures).   |            | We          | also      | searched       | the        |
Comparisonofrunningtimesofourtemplate-basedmethodwith STRING 49 database, which gives known and predicted
classicaldocking(usingPatchDock)ontheDockingBenchmark.(a)
|     |     | interactions | along | with | a   | confidence |     | score. | Here, all |
| --- | --- | ------------ | ----- | ---- | --- | ---------- | --- | ------ | --------- |
Runningtimesareplottedasafunctionofthenumberoftargetchains
|     |     | active prediction |     | methods |     | (neighborhood, |     | coexpression, |     |
| --- | --- | ----------------- | --- | ------- | --- | -------------- | --- | ------------- | --- |
inDockingBenchmark.Theanalysisisperformedontwotemplate
datasets(composedof1036and7922interfaces).(b)Runningtimeof gene fusion, co-occurence, coexpression, experiments,
ourtemplate-basedmethodisplottedasafunctionofthenumberof databases, homology, and text mining) are used to obtain
templateinterfacesfor165targetchainsinDockingBenchmark.If
|     |     | the confidence |     | score. | In  | this | way, | 34 interactions |     |
| --- | --- | -------------- | --- | ------ | --- | ---- | ---- | --------------- | --- |
therewerearound25,000templateinterfaces,twomethodswouldhave
thesamerunningtimesfor165targetproteins.[Colorfigurecanbe are found in STRING, besides those 30 interactions
viewedintheonlineissue,whichisavailableatwileyonlinelibrary.com.] (Table I). In this network, transcription factors such as
TableI
TheNumberofPredictedInteractionsandVerificationontheExperimentalData
|     | No.ofpredicted | No.ofverified |     | No.offurther |     |     |     | Totalno.ofverified |     |
| --- | -------------- | ------------- | --- | ------------ | --- | --- | --- | ------------------ | --- |
Templatedataset interactionsa interactionsb interactionsc interactions
Defaultd
| P53templates(default) | 108(49) | 30  |     |     | 34  |     |     |     | 64  |
| --------------------- | ------- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1036templates         | 53(38)  | 18  |     |     | 17  |     |     |     | 35  |
| Total                 | 114(50) | 31  |     |     | 37  |     |     |     | 68  |
Relaxede
| p53templates  | 396(68) | 52  |     |     | 105 |     |     |     | 157 |
| ------------- | ------- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1036templates | 721(71) | 67  |     |     | 177 |     |     |     | 244 |
| Total         | 870(71) | 83  |     |     | 222 |     |     |     | 305 |
aNumbersinparanthesesrepresentthenumberofproteins;thatis,108interactionsbetween49proteins.
bExperimentalinteractiondataforverificationareobtainedfromthehumanmolecularinteractionmap,DIP,MINT,BIND,andIntAct.
_ 49
cForfurtherevidenceforthepredictedinteractions,weusedtheSTRI NG searchtoolwhereweconsideredallactivepredictionmethods(neighborhood,coexpression,
genefusion,co-occurence,experiments,databases,homology,andtextmining)withmediumconfidencethreshold(0.4).
dThedefaultstructuralmatchingthresholdsareasfollows:40%oftheresiduesoftemplatechainsshouldgeometrically matchthetargetsurfacestopasstothenext
step.Thisthresholdis60%fortemplatechainscontaininglessthan50residues.
eIntherelaxedcase,defaultmatchingthresholdsarereducedby10%wherenewthresholdsare30and50%.
|     |     |     |     |     |     |     |     | PROTEINS | 1247 |
| --- | --- | --- | --- | --- | --- | --- | --- | -------- | ---- |

 10970134, 2012, 4, Downloaded from https://onlinelibrary.wiley.com/doi/10.1002/prot.24022 by Koc University, Wiley Online Library on [11/02/2026]. See the Terms and Conditions (https://onlinelibrary.wiley.com/terms-and-conditions) on Wiley Online Library for rules of use; OA articles are governed by the applicable Creative Commons License
N.Tuncbagetal.
E2F1-2-3, Max, Myc, Jun, and Fos interconnect via mul- similar interface architectures. Also, an accurate model
tiple interactions. As expected, there is also a large num- for a protein complex can be found from multiple
ber of interactions between cyclins and kinases. The tem- template interfaces. For instance, the matripase/trypsin
plate set containing 1036 interfaces gives just 53 interac- inhibitor is modeled on the template interfaces extracted
tions between 38 proteins with default thresholds, of from the thrombin/trypsin inhibitor complex, trypsino-
which 18 interactions have experimental evidence and 17 gen/trypsin inhibitor complex, and chymotrypsin/trypsin
49
| more are | verified | in STRING |     | (35 | verified | interactions |     | inhibitor | complex. |     |     |     |     |     |
| -------- | -------- | --------- | --- | --- | -------- | ------------ | --- | --------- | -------- | --- | --- | --- | --- | --- |
in total). When the matching thresholds are relaxed by In the validation, we used only rigid cases in the
10% (with new thresholds being 30 and 50%, respec- benchmark. As a limitation of our current method, it
tively), the template set in the MIM gives 396 putative does not perform well on difficult cases. Because the first
protein complexes between 68 protein chains of which a phase of this algorithm is rigid-body alignment, it is
total of 157 interactions are verified by interaction data- hard to handle conformational changes in the protein
bases. Using the relaxed matching thresholds, 721 inter- from the bound to the unbound states, if the conforma-
actions between 71 proteins are found from the 1036 tional change is in the binding region. As a future direc-
template interfaces of which a total of 244 interactions tion, we are working on integrating flexible aligment into
are verified (Table I). The results show that using strict the first phase of the method instead of rigid-body align-
matching thresholds not only give more reliable predic- ment and also integrating experimentally known multiple
tions but also miss true positives. When thresholds are conformations of template interfaces in the PDB. An
relaxed, the true positive rate increases; however, false additional key feature of this strategy is that the method
positives also increase. effectively distinguishes the nonbinders from binders.
|     |     |     |     |     |     |     |     | The verification |         | of the predicted | interactions |         | on an        | inde- |
| --- | --- | --- | --- | --- | --- | --- | --- | ---------------- | ------- | ---------------- | ------------ | ------- | ------------ | ----- |
|     |     |     |     |     |     |     |     | pendent          | protein | set (obtained    | from         | MIM)    | shows        | that  |
|     |     |     |     |     |     |     |     | this method      | can     | be used          | to predict   | whether | two proteins |       |
CONCLUSIONS
|     |     |     |     |     |     |     |     | interact, | in addition | to the | 3D modeling |     | of the | protein |
| --- | --- | --- | --- | --- | --- | --- | --- | --------- | ----------- | ------ | ----------- | --- | ------ | ------- |
Here, we presented the validation of a combinatorial complexes. Comparison of the running times with
approach to effectively model the 3D structures of pro- docking illustrates that the template-based approach is
tein complexes and the verification of the predicted pair- dramatically faster than docking. The limitation of this
wise interactions on a pathway scale. This approach relies template-based flexible docking approach is the diversity
on the expectation that the number of protein–protein of the template dataset. Currently, most of the available
interface architectures in nature is limited; thus, extrapo- high-resolution structures are of monomers, and the
lation of the known architecture space to target protein number of experimentally determined different interface
surfaces may help to identify protein interactions. This architectures is still limited. However, the projected fast
knowledge-based approach is made more physical by growth in the number of experimentally determined
combining it with flexible refinement of the solutions protein complexes in the near future will lead to an
and energy assessments to rank them. The docking increasing number of different interface architectures,
performance of the method is examined on the Docking which we expect to result in an increased use of such fast
Benchmark proteins. The validation results show that if a and reliable approaches by the community.
| structurally | similar | interface |     | is available |     | in the | template |     |     |     |     |     |     |     |
| ------------ | ------- | --------- | --- | ------------ | --- | ------ | -------- | --- | --- | --- | --- | --- | --- | --- |
dataset, the method can find the binding surface accu- ACKNOWLEDGMENTS
| rately and   | efficiently. |               | After self   | hits       | are        | eliminated  | from     |                   |              |            |                         |     |               |         |
| ------------ | ------------ | ------------- | ------------ | ---------- | ---------- | ----------- | -------- | ----------------- | ------------ | ---------- | ----------------------- | --- | ------------- | ------- |
|              |              |               |              |            |            |             |          | The authors       | thank        | Dr.        | Dina Duhovny-Schneidman |     |               | for     |
| the template |              | dataset,      | the method   |            | produces   | near        | native   |                   |              |            |                         |     |               |         |
|              |              |               |              |            |            |             |          | suggestions.      | The          | content    | of this publication     |     | necessarily   |         |
| models       | for 25       | complexes     | and          | acceptable |            | models      | for two  |                   |              |            |                         |     |               |         |
|              |              |               |              |            |            |             |          | neither           | does reflect | the views  | or policies             |     | of the        | Depart- |
| complexes    | out          | of 88         | with default |            | matching   | thresholds. |          |                   |              |            |                         |     |               |         |
|              |              |               |              |            |            |             |          | ment of           | Health       | and Human  | Services                | nor | does mention  |         |
| When the     | hot          | spot matching |              | threshold  |            | is relaxed, | this     |                   |              |            |                         |     |               |         |
|              |              |               |              |            |            |             |          | of trade          | names,       | commercial | products,               | or  | organizations |         |
| number       | increases    | to            | 41 near      | native     | models     |             | and 25   |                   |              |            |                         |     |               |         |
|              |              |               |              |            |            |             |          | imply endorsement |              | by the     | U.S. Government.        |     |               |         |
| acceptable   | models;      | however,      |              | the number |            | of false    | positive |                   |              |            |                         |     |               |         |
| models       | increases.   | Even          | if the       | sequence   | similarity |             | between  |                   |              |            |                         |     |               |         |
REFERENCES
| template   | interfaces | and    | target      | proteins | are         | decreased | dra-     |         |             |              |             |     |            |         |
| ---------- | ---------- | ------ | ----------- | -------- | ----------- | --------- | -------- | ------- | ----------- | ------------ | ----------- | --- | ---------- | ------- |
| matically, | the        | method | still finds | near     | native      | and       | accepta- |         |             |              |             |     |            |         |
|            |            |        |             |          |             |           |          | 1. Aloy | P, Bottcher | B, Ceulemans | H, Leutwein | C,  | Mellwig C, | Fischer |
| ble models | (17        | out of | 88 with     | default  | thresholds, |           | 61 out   |         |             |              |             |     |            |         |
S,GavinAC,BorkP,Superti-FurgaG,SerranoL,RussellRB.Struc-
| of 88 with | the     | hot spot   | threshold |     | relaxed | case). | The case |                     |          |            |           |     |           |         |
| ---------- | ------- | ---------- | --------- | --- | ------- | ------ | -------- | ------------------- | -------- | ---------- | --------- | --- | --------- | ------- |
|            |         |            |           |     |         |        |          | ture-based          | assembly | of protein | complexes |     | in yeast. | Science |
| studies    | provide | a detailed | view      | how | the     | method | predicts | 2004;303:2026–2029. |          |            |           |     |           |         |
accurate models. For example, the ribonuclease/barstar 2. KielC,BeltraoP,SerranoL.Analyzingproteininteractionnetworks
complex is modeled on a catalytic antibody 4B2 homo- usingstructuralinformation.AnnuRevBiochem2008;77:415–441.
|                 |     |                 |     |           |           |              |          | 3. Tuncbag    | N, Gursoy | A, Keskin | O. Prediction |     | of protein–protein  |     |
| --------------- | --- | --------------- | --- | --------- | --------- | ------------ | -------- | ------------- | --------- | --------- | ------------- | --- | ------------------- | --- |
| dimer. Although |     | the overall     |     | folds and | sequences |              | of these |               |           |           |               |     |                     |     |
|                 |     |                 |     |           |           |              |          | interactions: | unifying  | evolution | and structure | at  | protein interfaces. |     |
| two complexes   |     | are dissimilar, |     | both      | contain   | structurally |          |               |           |           |               |     |                     |     |
PhysBiol2011;8:035006.
1248 PROTEINS

Template-BasedDockingandFlexibleRefinement
4. Andrusier N, Mashiach E, Nussinov R, Wolfson HJ. Principles of 28. Aytuna AS, Gursoy A, Keskin O. Prediction of protein–protein
flexibleprotein–proteindocking.Proteins2008;73:271–289. interactions by combining structure and sequence conservation in
5. GrayJJ.High-resolutionprotein–proteindocking.CurrOpinStruct proteininterfaces.Bioinformatics2005;21):2850–2855.
Biol2006;16:183–193. 29. Ogmen U, Keskin O, Aytuna AS, Nussinov R, Gursoy A. PRISM:
6. Halperin I, Ma B, Wolfson H, Nussinov R. Principles of docking: protein interactions by structural matching. Nucleic Acids Res
an overview of search algorithms and a guide to scoring functions. 2005;33:W331–W336.
Proteins2002;47:409–443. 30. Gunther S, May P, Hoppe A, Frommel C, Preissner R. Docking
7. de Vries SJ, van Dijk M, Bonvin AM. The HADDOCKweb server without docking: ISEARCH—prediction of interactions using
fordata-drivenbiomoleculardocking.NatProtoc2010;5:883–897. knowninterfaces.Proteins2007;69:839–844.
8. LeskVI,SternbergMJ.3D-Garden:asystemformodellingprotein– 31. Sinha R, Kundrotas PJ, Vakser IA. Docking by structural similarity
protein complexes based on conformational refinement of ensem- atprotein–proteininterfaces.Proteins2010;78:3235–3241.
bles generated with the marching cubes algorithm. Bioinformatics 32. Aramini JM, Ma LC, Zhou L, Schauder CM, Hamilton K, Amer
2008;24:1137–1144. BR, Mack TR, Lee HW, Ciccosanti CT, Zhao L, Xiao R, Krug RM,
9. Cheng TM, Blundell TL, Fernandez-Recio J. pyDock: electrostatics Montelione GT. Dimer interface of the effector domain of non-
and desolvation for effective scoring of rigid-body protein–protein structuralprotein1frominfluenzaAvirus:aninterfacewithmulti-
docking.Proteins2007;68:503–515. plefunctions.JBiolChem2011;286:26050–26060.
10. Janin J. Protein–protein docking tested in blind predictions: the 33. Kundrotas PJ, Vakser IA. Accuracy of protein–protein binding sites
CAPRIexperiment.MolBiosyst2010;6:2351–2362. in high-throughput template-based modeling. PLoS Comput Biol
11. Wodak SJ, Mendez R. Prediction of protein–protein interactions: 2010;6:e1000727.
the CAPRI experiment, its evaluation and implications. Curr Opin 34. Tuncbag N, Gursoy A, Nussinov R, Keskin O. Predicting protein–
StructBiol2004;14:242–249. protein interactions on a proteome scale by matching evolutionary
12. Kastritis PL, Bonvin AM. Are scoring functions in protein–protein and structural similarities at interfaces using PRISM. Nat Protoc
docking ready to predict interactomes? Clues from a novel binding 2011;6:1341–1354.
affinitybenchmark.JProteomeRes2010;9:2216–2225. 35. Mashiach E, Nussinov R, Wolfson HJ. FiberDock: flexible induced-
13. FeliuE,AloyP,OlivaB.Ontheanalysisofprotein–proteininterac- fit backbone refinement in molecular docking. Proteins 2010;78:
tions via knowledge-based potentials for the prediction of protein– 1503–1519.
proteindocking.ProteinSci2011;20:529–541. 36. Kohn KW. Molecular interaction map of the mammalian cell
14. Tsai CJ, Lin SL, Wolfson HJ, Nussinov R. A dataset of protein– cyclecontrolandDNArepairsystems.MolBiolCell1999;10:2703–
protein interfaces generated with a sequence-order-independent 2734.
comparisontechnique.JMolBiol1996;260:604–620. 37. Hwang H, Pierce B, Mintseris J, Janin J, Weng Z. Protein–protein
15. TsaiCJ,LinSL,WolfsonHJ,NussinovR.Protein–proteininterfaces: dockingbenchmarkversion3.0.Proteins2008;73:705–709.
architectures and interactions in protein–protein interfaces and in 38. TuncbagN,KeskinO,GursoyA.HotPoint:hotspotpredictionserver
protein cores. Their similarities and differences. Crit Rev Biochem forproteininterfaces.NucleicAcidsRes2010;38(suppl):W402–W406.
MolBiol1996;31:127–152. 39. Hubbard SJ, Thornton JM. NACCESS. University College, London:
16. Tsai CJ, Lin SL, Wolfson HJ, Nussinov R. Studies of protein– DepartmentofBiochemistryandMolecularBiology; 1993.
protein interfaces: a statistical analysis of the hydrophobic effect. 40. Nussinov R, Wolfson HJ. Efficient detection of three-dimensional
ProteinSci1997;6:53–64. structural motifs in biological macromolecules by computer vision
17. Tsai CJ, Xu D, Nussinov R. Structural motifs at protein–protein techniques.ProcNatlAcadSciUSA1991;88:10495–10499.
interfaces: protein cores versus two-state and three-state model 41. Shatsky M, Nussinov R, Wolfson HJ. A method for simultaneous
complexes.ProteinSci1997;6:1793–1805. alignmentofmultipleproteinstructures.Proteins2004;56:143–156.
18. Tsai CJ, Xu D, Nussinov R. Protein folding via binding and vice 42. Chailyan A, Marcatili P, Tramontano A. The association of heavy
versa.FoldDes1998;3:R71–R80. and light chain variable domains in antibodies: implications for
19. Tuncbag N, Gursoy A, Guney E, Nussinov R, Keskin O. Architec- antigenspecificity.FEBSJ2011;278:2858–2866.
tures and functional coverage of protein–protein interfaces. J Mol 43. Gao M, Skolnick J. New benchmark metrics for protein–protein
Biol2008;381:785–802. dockingmethods.Proteins2011;79:1623–1634.
20. Keskin O, Ma B, Rogale K, Gunasekaran K, Nussinov R. Protein– 44. Schneidman-Duhovny D, Inbar Y, Nussinov R, Wolfson HJ. Patch-
protein interactions: organization, cooperativity and mapping in a Dock and SymmDock: servers for rigid and symmetric docking.
bottom-upsystemsbiologyapproach.PhysBiol2005;2:S24–S35. NucleicAcidsRes2005;33:W363–W367.
21. Keskin O, Nussinov R. Favorable scaffolds: proteins with different 45. SalwinskiL,MillerCS,SmithAJ,PettitFK,BowieJU,EisenbergD.
sequence, structure and function may associate in similar ways. The database of interacting proteins: 2004 update. Nucleic Acids
ProteinEngDesSel2005;18:11–24. Res2004;32:D449–D451.
22. AloyP,RussellRB.Interrogatingproteininteractionnetworksthrough 46. Chatr-aryamontriA,CeolA,PalazziLM,NardelliG,SchneiderMV,
structuralbiology.ProcNatlAcadSciUSA2002;99:5896–5901. Castagnoli L, Cesareni G. MINT: the Molecular INTeraction data-
23. Bell RE, Ben-Tal N. In silico identification of functional protein base.NucleicAcidsRes2007;35:D572–D574.
interfaces.CompFunctGenomics2003;4:420–423. 47. BaderGD,BetelD,HogueCW.BIND:thebiomolecularinteraction
24. GlaserF,PupkoT,PazI,BellRE,Bechor-ShentalD,MartzE,Ben-TalN. networkdatabase.NucleicAcidsRes2003;31:248–250.
ConSurf:identificationoffunctionalregionsinproteinsbysurface-map- 48. Kerrien S, Alam-Faruque Y, Aranda B, Bancarz I, Bridge A, Derow
pingofphylogeneticinformation.Bioinformatics2003;19:163–164. C,DimmerE,FeuermannM,FriedrichsenA,HuntleyR,KohlerC,
25. KeskinO,MaB,NussinovR.Hotregionsinprotein–proteininter- Khadake J, Leroy C, Liban A, Lieftink C, Montecchi-Palazzi L,
actions:theorganizationandcontributionofstructurallyconserved OrchardS,RisseJ,RobbeK,RoechertB,ThorneycroftD,ZhangY,
hotspotresidues.JMolBiol2005;345:1281–1294. Apweiler R, Hermjakob H. IntAct—open source resource for
26. KeskinO,TsaiCJ,WolfsonH,NussinovR.Anew,structurallynon- molecularinteractiondata.NucleicAcidsRes2007;35:D561–D565.
redundant, diverse data set of protein–protein interfaces and its 49. SzklarczykD,FranceschiniA,KuhnM,SimonovicM,RothA,Min-
implications.ProteinSci2004;13:1043–1055. guezP,DoerksT,StarkM,MullerJ,BorkP,JensenLJ,vonMering
27. Keskin O, Gursoy A, Ma B, Nussinov R. Principles of protein– C. The STRING database in 2011: functional interaction networks
protein interactions: what are the preferred ways for proteins to of proteins, globally integrated and scored. Nucleic Acids Res
interact?ChemRev2008;108:1225–1244. 2011;39:D561–D568.
PROTEINS 1249
10970134,
2012,
4,
Downloaded
from
https://onlinelibrary.wiley.com/doi/10.1002/prot.24022
by
Koc
University,
Wiley
Online
Library
on
[11/02/2026].
See
the
Terms
and
Conditions
(https://onlinelibrary.wiley.com/terms-and-conditions)
on
Wiley
Online
Library
for
rules
of
use;
OA
articles
are
governed
by
the
applicable
Creative
Commons
License