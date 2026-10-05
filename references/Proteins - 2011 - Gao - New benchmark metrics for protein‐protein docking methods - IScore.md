proteins
STRUCTUREOFUNCTIONOBIOINFORMATICS
| New     | benchmark   |           |     |     | metrics |     | for |     | protein-protein |     |     |     |     |     |     |
| ------- | ----------- | --------- | --- | --- | ------- | --- | --- | --- | --------------- | --- | --- | --- | --- | --- | --- |
| docking |             | methods   |     |     |         |     |     |     |                 |     |     |     |     |     |     |
| Mu Gao  | and Jeffrey | Skolnick* |     |     |         |     |     |     |                 |     |     |     |     |     |     |
CenterfortheStudyofSystemsBiology,SchoolofBiology,GeorgiaInstituteofTechnology,Atlanta,Georgia30318
INTRODUCTION
ABSTRACT
|     |     |     |     |     | In the | quest to | determine | all | protein-protein |     | interactions |     | in  | a given | pro- |
| --- | --- | --- | --- | --- | ------ | -------- | --------- | --- | --------------- | --- | ------------ | --- | --- | ------- | ---- |
With the development of many compu- teome, recent high-throughput technologies have enabled substantial pro-
| tational     | methods that       | predict  | the struc- |         | 1–3 |        |               |         |         |         |                 |           |        |               |      |
| ------------ | ------------------ | -------- | ---------- | ------- | --- | ------ | ------------- | ------- | ------- | ------- | --------------- | --------- | ------ | ------------- | ---- |
|              |                    |          |            | gress.  |     | Drafts | for different |         | model   | systems | are             | emerging, | though |               | many |
| tural models | of protein-protein |          | com-       |         |     |        |               |         |         |         |                 |           |        |               |      |
|              |                    |          |            | details | are | still  | missing.      | 4–7 The | mapping | of      | protein-protein |           |        | interactions, |      |
| plexes,      | there is a         | pressing | need to    |         |     |        |               |         |         |         |                 |           |        |               |      |
benchmark their performance. As was however, is just a starting point toward revealing their functional roles in liv-
|          |             |           |            | ing    | biosystems. | In           | order        | to understand |                    | protein-protein |         | interactions, |           | it is | nec- |
| -------- | ----------- | --------- | ---------- | ------ | ----------- | ------------ | ------------ | ------------- | ------------------ | --------------- | ------- | ------------- | --------- | ----- | ---- |
| the case | for protein | monomers, | assess-    |        |             |              |              |               |                    |                 |         |               |           |       |      |
|          |             |           |            | essary | to          | structurally | characterize |               | all representative |                 | protein |               | complexes | at    | high |
| ing the  | quality of  | models    | of protein |        |             |              |              |               |                    |                 |         |               |           |       |      |
8
| complexes | is not straightforward. |     | An  | resolution. |     |     |     |     |     |     |     |     |     |     |     |
| --------- | ----------------------- | --- | --- | ----------- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
effective scoring scheme should be able Despite rapid growth in the number of structurally solved protein com-
| to detect | substructure    | similarity    | and   |           | 9   |                    |              |               |               |     |           |        |      |          |     |
| --------- | --------------- | ------------- | ----- | --------- | --- | ------------------ | ------------ | ------------- | ------------- | --- | --------- | ------ | ---- | -------- | --- |
|           |                 |               |       | plexes,   |     | the pace           | of structure |               | determination |     | lags far  | behind | the  | pace of  | the |
| estimate  | its statistical | significance. | Here, |           |     |                    |              |               |               |     |           |        |      |          |     |
|           |                 |               |       | detection |     | of protein-protein |              | interactions. |               | To  | fill this | gap,   | many | computa- |     |
wefocusoncharacterizingthesimilarity tional approaches have been proposed for predicting the structures of protein
| of the interfaces  | of          | the complex | and         |            |     |               |        |             |             |                  |            |        |                |         |     |
| ------------------ | ----------- | ----------- | ----------- | ---------- | --- | ------------- | ------ | ----------- | ----------- | ---------------- | ---------- | ------ | -------------- | ------- | --- |
|                    |             |             |             | complexes. |     | They          | can be | roughly     | categorized | into             | two        | types: | Template-Based |         |     |
| introduce          | two scoring | functions.  | The         |            |     |               |        |             |             |                  |            |        |                |         |     |
|                    |             |             |             | (TB)       | and | Template-Free |        | (TF).       | In TB       | approaches,10–15 |            |        | one first      | builds  | a   |
| first, the         | interfacial | Template    | Modeling    |            |     |               |        |             |             |                  |            |        |                |         |     |
|                    |             |             |             | homology   |     | model         | based  | on a solved | template    |                  | structure, | and    | then           | refines | the |
| score (iTM-score), | measures    |             | the geomet- |            |     |               |        |             |             |                  |            |        |                |         |     |
16–23
|              |             |           |             | model.   |     | In TF approaches, |         |     | also       | known | as        | protein-protein |        | docking  |     |
| ------------ | ----------- | --------- | ----------- | -------- | --- | ----------------- | ------- | --- | ---------- | ----- | --------- | --------------- | ------ | -------- | --- |
| ric distance | between     | the       | interfaces, |          |     |                   |         |     |            |       |           |                 |        |          |     |
|              |             |           |             | methods, |     | one docks         | unbound |     | components |       | that form | the             | target | complex. |     |
| while the    | second, the | Interface | Similar-    |          |     |                   |         |     |            |       |           |                 |        |          |     |
ity score (IS-score), evaluates their Both methods have advantages and disadvantages. TB approaches generally
residue-residue contact similarity in have higher accuracy, but suffer from low coverage because of their depend-
addition to their geometric similarity. ence on the availability of template structures. Although the issue of low
We first demonstrate that the IS-score is coverage might be overcome by the recognition that the structural space of
more suitable for assessing docking protein-protein interfaces is highly degenerate,24 in practice identifying
models thantheiTM-score.TheIS-score
|                   |             |               |            | which   | protein        | pairs      | actually | interact     | is           | very challenging. |             | On    | the             | other    | hand, |
| ----------------- | ----------- | ------------- | ---------- | ------- | -------------- | ---------- | -------- | ------------ | ------------ | ----------------- | ----------- | ----- | --------------- | -------- | ----- |
| is then validated | in          | a large-scale | bench-     |         |                |            |          |              |              |                   |             |       |                 |          |       |
|                   |             |               |            | TF      | approaches     | can        | deal     | with a       | novel target | whose             | quaternary  |       | structure       |          | does  |
| mark test         | on 1562     | dimeric       | complexes. |         |                |            |          |              |              |                   |             |       |                 |          |       |
|                   |             |               |            | not     | match          | any solved | template |              | structure,   | but               | there       | is no | guarantee       | of       | high- |
| Finally,          | the scoring | function      | is applied |         |                |            |          |              |              |                   |             |       |                 |          |       |
|                   |             |               |            | quality | docking        | models,    |          | particularly | when         | bound             | structures  |       | undergo         | signifi- |       |
| to evaluate       | docking     | models        | submitted  |         |                |            |          |              |              |                   |             |       |                 |          |       |
|                   |             |               |            | cant    | conformational |            | changes  | from         | the          | unbound           | structures. |       | 25 Furthermore, |          |       |
| to the Critical   | Assessment  | of            | Prediction |         |                |            |          |              |              |                   |             |       |                 |          |       |
of Interactions (CAPRI) experiments. TF approaches require the information that two input proteins interact; that
While the results according to the new is, they are not reliable in predicting whether two proteins interact or not,
|                |     |           |            | largely | due | to the | limitations |     | of force | fields | used for | evaluating |     | interaction |     |
| -------------- | --- | --------- | ---------- | ------- | --- | ------ | ----------- | --- | -------- | ------ | -------- | ---------- | --- | ----------- | --- |
| scoring scheme | are | generally | consistent |         |     |        |             |     |          |        |          |            |     |             |     |
energy.26
with the original CAPRI assessment, the By comparison, TB methods usually contain (explicitly or implic-
IS-score identifies models whose signifi- itly) an evolutionary component, which prefers templates sharing conserved
| cancewas | previously | underestimated. |     |            |           |              |         |               |           |            |     |          |      |            |     |
| -------- | ---------- | --------------- | --- | ---------- | --------- | ------------ | ------- | ------------- | --------- | ---------- | --- | -------- | ---- | ---------- | --- |
|          |            |                 |     | biological |           | interactions | with    | target        | proteins. | Thus,      | in  | addition | to   | predicting |     |
|          |            |                 |     | the        | structure | of           | protein | interactions, |           | TB methods |     | may be   | used | to predict |     |
Proteins2011;79:1623–1634.
| VVC 2011Wiley-Liss,Inc. |     |     |     | whether |              | two proteins | interact.       |     |     |         |          |     |                |     |     |
| ----------------------- | --- | --- | --- | ------- | ------------ | ------------ | --------------- | --- | --- | ------- | -------- | --- | -------------- | --- | --- |
|                         |     |     |     |         | To benchmark |              | the performance |     | of  | docking | methods, | a   | community-wide |     |     |
Keywords: docking; protein-protein experiment, known as CAPRI, has been carried out.27–29 One central task is
| interaction; | protein-protein |           | interface; |     |     |     |     |     |     |     |     |     |     |     |     |
| ------------ | --------------- | --------- | ---------- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| structure    | prediction;     | TM-score; | IS-score;  |     |     |     |     |     |     |     |     |     |     |     |     |
CAPRI. Grantsponsor:NationalInstitutesofHealth;Grantnumber:GM-48835.
|     |     |     |     | *Correspondence |     | to: Jeffrey | Skolnick, | Center | for the | Study of Systems | Biology, | School | of  | Biology, | Georgia |
| --- | --- | --- | --- | --------------- | --- | ----------- | --------- | ------ | ------- | ---------------- | -------- | ------ | --- | -------- | ------- |
InstituteofTechnology,25014thStreetNW,Atlanta,GA30318.E-mail:skolnick@gatech.edu.
Received8November2010;Revised22December2010;Accepted30December2010
Publishedonline18January2011inWileyOnlineLibrary(wileyonlinelibrary.com).DOI:10.1002/prot.22987
| VVC 2011WILEY-LISS,INC. |     |     |     |     |     |     |     |     |     |     |     |     | PROTEINS |     | 1623 |
| ----------------------- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | -------- | --- | ---- |

 10970134, 2011, 5, Downloaded from https://onlinelibrary.wiley.com/doi/10.1002/prot.22987 by Koc University, Wiley Online Library on [16/02/2026]. See the Terms and Conditions (https://onlinelibrary.wiley.com/terms-and-conditions) on Wiley Online Library for rules of use; OA articles are governed by the applicable Creative Commons License
M.GaoandJ.Skolnick
to measure the quality of a predicted docking model, the IS-score to docking models submitted to the CAPRI
| using its      | target structure,  | usually | a            | solved   | crystal     | struc- experiments. |     |     |     |
| -------------- | ------------------ | ------- | ------------ | -------- | ----------- | ------------------- | --- | --- | --- |
| ture, as       | the gold standard. |         | Furthermore, |          | in the case | of                  |     |     |     |
| template-based | modeling,          | it      | is also      | critical | to measure  |                     |     |     |     |
METHODS
| the quality | of both | the template | and | the | final | model. |     |     |     |
| ----------- | ------- | ------------ | --- | --- | ----- | ------ | --- | --- | --- |
A˚
Thus, any improvement or deterioration resulting from A heavy-atom distance cutoff of 4.5 is employed to
the ‘‘refinement’’ procedure, designed to improve over define an interfacial contact. A protein-protein interface
| the template | alignment, | can | be evaluated. |     | For these | pur-              |                 |               |              |
| ------------ | ---------- | --- | ------------- | --- | --------- | ----------------- | --------------- | ------------- | ------------ |
|              |            |     |               |     |           | is the collection | of all residues | with at least | one interfa- |
poses, one needs to derive effective structure comparison cial contact between pairs of proteins.
| metrics. | The CAPRI | assessors | employed | complex | criteria |     |     |     |     |
| -------- | --------- | --------- | -------- | ------- | -------- | --- | --- | --- | --- |
based on the Root Mean Square Deviation (RMSD) and Scoringfunctionandsearchalgorithm
29
| the fraction | of conserved |     | native contacts |     | f . | While |     |     |     |
| ------------ | ------------ | --- | --------------- | --- | --- | ----- | --- | --- | --- |
nat
these criteria are convenient, they have three limitations: Assuming that a native (target) structure has L interfa-
|             |                |          |              |              |           | cial residues,         | the iTM-score        | of a corresponding | docking       |
| ----------- | -------------- | -------- | ------------ | ------------ | --------- | ---------------------- | -------------------- | ------------------ | ------------- |
| The first   | is that RMSD   | is often | dominated    |              | by the    | largest                |                      |                    |               |
|             |                |          |              |              |           | model is               | defined by comparing | the geometric      | distances     |
| deviations, | and hence,     | may      | overlook     | substructure | similar-  |                        |                      |                    |               |
|             |                |          |              |              |           | of the native          | interfacial          | residues of the    | model and the |
| ity. The    | second is that | the      | statistical  | significance |           | of a                   |                      |                    |               |
|             |                |          | dependent.30 |              |           | native structure,33,36 |                      |                    |               |
| given RMSD  | value is       | length   |              |              | The third | is                     |                      |                    |               |
that the thresholds employed for model quality classifica- " #
XNa
1
tion are often subjective, in the sense that an assessment iTM(cid:1)score¼ max 1=ð1þd2=d2Þ ð1Þ
i 0
of the statistical significance of the given structural com- L
i¼1
| parison | metric is lacking. |     |     |     |     |     |     |     |     |
| ------- | ------------------ | --- | --- | --- | --- | --- | --- | --- | --- |
The problem of model quality assessment is not where N is the number of superimposed native interfa-
a
unique to protein docking experiments. An analogous cial residues, d is the Euclidean distance between the Ca
i
situation was encountered in the evaluation of structural atoms from the ith superimposed repsidffiffiuffiffiffieffiffiffiffipffiffiffiaffiffiiffirffi, and the
models predicted for monomeric proteins. In the recent empirical scaling factor d (cid:2)1:243 L (cid:1)15(cid:1)1:8 is
|     |     |     |     |     |     |     |     | 0   | Q   |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
Critical Assessment of Protein Structure Prediction introduced to correct for length effects. Note that the
(CASP), several commonly used scoring functions definition of the iTM-score is exactly the same as used
include the Global Distance Test (GDT) score, 31 the for assessing the model quality of the global structure
|     | 32  |     |     | 33  |     |     |     | 33  |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
MaxSub score, and the TM-score. The statistical sig- alignment of monomeric proteins. However, the TM-
nificance relative to random of both the GDT and the scores of interfaces and of individual proteins have a dif-
MaxSub scores are sensitive to the size of the target pro- ferent level of statistical significance at the same numeri-
33
tein. As a result, one often cannot tell whether a raw cal value (see below). To avoid confusion, we use the
score indicates a significant prediction. By contrast, the term iTM-score to denote the TM-score of interfaces and
TM-score corrects for length effects. Based on the statis- reserve the notation TM-score for the global comparison
tics obtained from comparing random protein structures of a pair of structures.
at various lengths, a TM-score of 0.4 or higher indicates In order to calculate the distance d, a subset of corre-
i
33
a significant prediction. Other statistically rigorous sponding residues are superimposed using the Kabsch
treatments also have been undertaken to calculate the algorithm, 37 which minimizes their pairwise root-mean-
34,35
significance (i.e., P values) of protein models. square deviation, RMSD. Since there are many ways to
Previously, we introduced the iTM-score and the IS- select the subset, the notation max in Eq. (1) indicates
scoreiniAlign,36aprogramforthestructuralcomparison that the iTM-score is the maximum out of all possible
of protein-protein interfaces based on interface structure superimpositions. A heuristic iterative extension algo-
alignments, where the equivalence of target and template rithm is employed to calculate the iTM-score,33 similar
31
residues is not a priori specified. It has been shown that to the one used for calculating the GDT-score and
|     |     |     |     |     |     | 32  |     |     | 5   |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
the IS-score is an effective metric for evaluating structural MaxSub. Briefly, we select fragments of size L L,
sub
alignments of protein-protein interfaces. 24,36 In this L/2, L/4, ..., 4, respectively. When L is less than L,
sub
study, we examine both scoring functions for measuring initial fragments are selected by sliding continuously
the qualityof docking models. The keydifferencebetween along the native interface. Starting from an initial frag-
the previous study and the current one is that iAlign does ment of size L sub , the corresponding residues within L sub
not require any previously specified sequence correspon- in the model and native interfaces are superimposed.
dence, whereas in the current scenario, the mapping of Then, all model/interface residue pairs within a distance
equivalent target-template residuesis specifiedin advance. less than d are collected and superimposed again. The
0
As a result, one needs to adjust the random background process is iterated until the rigid-body transformation
| and recalibrate | the statistical |     | models, | as  | detailed | below. converges. |     |     |     |
| --------------- | --------------- | --- | ------- | --- | -------- | ----------------- | --- | --- | --- |
Furthermore, we performed large-scale benchmark tests The second scoring function is the Interface Similarity
to compare and validate our scoring schemes and applied score (IS-score), which measures not only geometric
1624 PROTEINS

 10970134, 2011, 5, Downloaded from https://onlinelibrary.wiley.com/doi/10.1002/prot.22987 by Koc University, Wiley Online Library on [16/02/2026]. See the Terms and Conditions (https://onlinelibrary.wiley.com/terms-and-conditions) on Wiley Online Library for rules of use; OA articles are governed by the applicable Creative Commons License
AssessmentofProteinDockingModels
Figure 1
Distributionsofrandomlyselectedprotein-proteininterfacepairs.(A)MeanofIS-scoresandinterfacialRMSDvaluesversusthesizeofprotein
interfaces.Horizontaldashedlinesarelocatedat0.1.(B)DistributionoftheIS-scoreamongrandominterfaces.Thehistogramistheobserved
scoredistribution,andthesolidlineisthefitaccordingtotheGumbeldistribution[Eq.(5)].[Colorfigurecanbeviewedintheonlineissue,
whichisavailableatwileyonlinelibrary.com.]
contacts.36
distances but also the conservation of interfacial Statisticalsignificance
| The IS-score | is derived | from the iTM-score | as follows, |     |                 |     |              |                 |              |
| ------------ | ---------- | ------------------ | ----------- | --- | --------------- | --- | ------------ | --------------- | ------------ |
|              |            |                    |             |     | The statistical |     | significance | of the IS-score | is estimated |
IS(cid:1)score¼ðSþs Þ=ð1þs Þ ð2Þ by comparing 24,120 randomly selected interface pairs of
|     |     | 0   | 0   |     |              |          |            |               |                   |
| --- | --- | --- | --- | --- | ------------ | -------- | ---------- | ------------- | ----------------- |
|     |     |     |     |     | same lengths | (see     | Data Set). | In each pair, | interfacial resi- |
|     |     | "   | #   |     | dues at      | the same | positions  | in respective | sequences are     |
XNa
|     | 1   |           |     |     | arbitrarily | assigned      | as equivalent. | Figure    | 1(A) shows the      |
| --- | --- | --------- | --- | --- | ----------- | ------------- | -------------- | --------- | ------------------- |
|     | ¼   | =ð1þd 2=d | 2Þ  | ð3Þ |             |               |                |           |                     |
|     | S   | max f     |     |     | means of    | the IS-scores |                | and iRMSD | values of unrelated |
|     | L   | i i       | 0   |     |             |               |                |           |                     |
i¼1
|     |     |     |     |     | interfaces. | Without    | applying | the scaling       | factor, the raw |
| --- | --- | --- | --- | --- | ----------- | ---------- | -------- | ----------------- | --------------- |
|     |     |     |     |     | IS-score    | calculated | using    | Eq. (3) decreases | exponentially   |
(cid:2)ðc=a þc=bÞ=2,
Here, the contact overlap factor f i i i i i as the length of the interface increases. Likewise, the
where a is the number of interfacial contacts observed at mean random iRMSD value increases exponentially as
i
the ith position of the native interface, b is the number the interface size increases. By comparison, the rescaled
i
of interfacial contacts observed at the corresponding IS-scores are approximately length-independent at a
position in the model, and c is the number of interfacial mean value of 0.10. It should be noted that the mean of
i
5
contacts conserved in both interfaces. If c i 0, f i is 0, random IS-scores calculated here is smaller than the
regardless of the value of b. The scaling factor mean of random IS-scores calculated previously with the
i
s (cid:2)0:14(cid:1)0:2=L0:3 is introduced to make the means of program iAlign. 36 The reason is that iAlign does not a
| 0   | Q   |     |     |     |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
the IS-scores length-independent among randomly priori impose a one-to-one sequence correspondence.
selected interfaces (see below). Note that the scaling fac- Therefore, iAlign usually finds a better correspondence
|                 |           |               |                    |     | (or alignment), |     | which gives | a higher | IS-score even for |
| --------------- | --------- | ------------- | ------------------ | --- | --------------- | --- | ----------- | -------- | ----------------- |
| tor is slightly | different | from what was | derived previously |     |                 |     |             |          |                   |
iAlign.36
in The adjustment is introduced to correct for randomly related interfaces.
a small shift in the means of the IS-scores among ran- Since the IS-scores are maxima, the extreme value dis-
dom interfaces. The search algorithm for calculating the tribution is a suitable statistical model for describing
IS-score is essentially the same as describe above for the their distribution. As shown in Figure 1(B), the probabil-
iTM-score. ity density function of the IS-scores calculated from
Both the iTM/IS-score give a maximum score of one the random background follows the extreme value
| for a perfect | model. |     |     |     | distribution, |     |     |     |     |
| ------------- | ------ | --- | --- | --- | ------------- | --- | --- | --- | --- |
PROTEINS 1625

M.GaoandJ.Skolnick
TableI complex are termed as the ligand/receptor of the com-
StatisticalSignificanceoftheIS-ScoresDerivedfrom24,120Pairsof plex. Let N denote the number of interfacial contacts
c
RandomInterfaces
observed in the native complex structure, and n the
IS-score number of native interfacial contacts preserved in the
docking model. The fraction of native contacts is f :
Pvalue Model Empirical nat
n/N. The interfacial RMSD, iRMSD, is the RMSD of the
c
0.05 0.125 0.120
0.01 0.142 0.134
Ca atoms of interfacial residues observed in a native
structure with respect to their positions in a docking
0.005 0.149 0.141
0.001 0.166 0.160 model, and the ligand RMSD, lRMSD, is the global
1e-04 0.190 0.187 RMSD of the Ca atoms of all ligand residues. The
1e-05 0.214 —
iRMSD is calculated after superimposing these native
1e-06 0.238 —
1e-08 0.286 — interfacial residues, whereas the lRMSD is calculated after
1e-10 0.334 — superimposing the receptors.
Datasets
Randombackground
fðzÞ¼exp½z(cid:1)expðzÞ(cid:3) ð4Þ
The random background for statistical significance
analysis was derived from the M-TASSER template
where z denotes the Z-score given by z 5 (s 2 l)/r. The library. 11 We first obtained all-against-all pairs of all di-
variable s denotes the IS-score; l is the location parame- meric complexes. A pair of dimers was then selected, if
ter, and r is the scale parameter. The corresponding P any two monomers, one from each dimer, have a global
value of the score can be calculated according to the sequence identity <30% and a global TM-score <0.4.
formula This selection led to a set of globally unrelated dimer
pairs. Since IS-score requires that the two interfaces have
P ¼1(cid:1)exp½(cid:1)expð(cid:1)zÞ(cid:3) ð5Þ
the same length, we randomly removed interfacial resi-
dues of the longer interface, if the two interfaces are of
The scores from random interfaces were fit to Eq. (4).
different size. The removal was carefully done by requir-
The resulting P values and their corresponding IS-scores
ing that all remaining interfacial residues maintain at
are given in Table I. The calculated P values according to
least one interfacial contact. To prevent possible over-rep-
the statistical model agree with the empirical values
resentation of any given dimer, we further required that
obtained by ranking the IS-scores of 24,120 random
no dimer appears more than 20 times in the final selec-
interface pairs. One may use these scores to quickly esti-
tions. The procedure yielded 24,120 pairs of interfaces,
mate statistical significance.
which were used for estimating the statistical significance
An improved estimation of statistical significance is
of the IS-score. In each pair, two interfacial residues were
obtained by modeling the distributions of scores at spe-
assigned as equivalent if they appear at the same posi-
cific lengths. Figure 2 shows the observed and modeled
tions in respective sequences after removing all non-
distributions at various lengths. Each distribution is
interfacial residues.
modeled by the Gumbel distribution described in Eq.
(4). The location and scale parameters can be estimated
Decoyset
through linear regression fits,
For the comparison between the iTM-score and the
l¼aþblnðLÞ 38 ð6Þ IS-score, we used a decoy set from the Dockgound.
r¼cþdlnðLÞ The decoy set was curated from docking models gener-
ated with unbound protein structures for 61 target com-
The parameters a to d, given in Table II, were obtained plexes. We further define a near native docking model if
by linear fitting to the location and scale parameters, it has lRMSD (cid:4)5 A˚ and f nat > 30%, and define an incor-
which were obtained through maximum likelihood esti- rect model if it has lRMSD >5 A˚ and f nat 5 0%. The
mates with the EVD package in the statistical platform R procedure produced 425 near native models and 5,232
(available at: http://www.r-project.org/). incorrect models.
Analysismeasures Dockingset
11
In addition to the iTM/IS-score, we also define com- From the M-TASSER template library, we selected
29 mon metrics adopted for evaluating docking models. 1,526 complexes whose individual proteins are less than
The smaller/larger of the two monomers in a binary 500 amino acids in length. Rigid-body docking using the
1626 PROTEINS
10970134,
2011,
5,
Downloaded
from
https://onlinelibrary.wiley.com/doi/10.1002/prot.22987
by
Koc
University,
Wiley
Online
Library
on
[16/02/2026].
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

 10970134, 2011, 5, Downloaded from https://onlinelibrary.wiley.com/doi/10.1002/prot.22987 by Koc University, Wiley Online Library on [16/02/2026]. See the Terms and Conditions (https://onlinelibrary.wiley.com/terms-and-conditions) on Wiley Online Library for rules of use; OA articles are governed by the applicable Creative Commons License
AssessmentofProteinDockingModels
Figure 2
DistributionsoftheZ-scoreamongrandominterfacesofvariouslengths.Longdashedlinesaretheobservedprobabilitydensity,andtheshort
dashedlinesaredirectfitsusingtheGumbeldistributions.SolidlinesareprobabilitydensitiescalculatedfortheIS-scoreswithstatisticalmodels
describedbyEqs.(4)and(6).L representsthelengthofquery,andN isthenumberofsamples.[Colorfigurecanbeviewedintheonlineissue,
|     |     | T   |     |     | S   |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
whichisavailableatwileyonlinelibrary.com.]
| bound structures       | from         | the complexes   | were    | subsequently |         | CAPRImodels                           |          |              |                 |      |
| ---------------------- | ------------ | --------------- | ------- | ------------ | ------- | ------------------------------------- | -------- | ------------ | --------------- | ---- |
| carried                | out with the | program FT-Dock | 23      | using        | default |                                       |          |              |                 |      |
|                        |              |                 |         |              |         | The docking                           | models   | for recent   | CAPRI targets   | were |
| parameters.            | The top      | 100 docking     | models, | ranked       | by      |                                       |          |              |                 |      |
|                        |              |                 |         |              |         | downloaded                            | from the | official web | site (available | at:  |
| shape complementarity, |              | were retained   | for     | validating   | the     |                                       |          |              |                 |      |
|                        |              |                 |         |              |         | http://www.ebi.ac.uk/msd-srv/capri/). |          |              | We selected     | ten  |
statistical significance of the IS-scores. In total, we col- recent protein-protein targets (T242T36, except for can-
lected 152,600 models by pooling together the top 100 celled T26, and RNA/protein targets T33 and T34), for
docking models from all complexes. which the docking models were available to the public.
|     |     |     |     |     |     |     |     |     | PROTEINS | 1627 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | -------- | ---- |

 10970134, 2011, 5, Downloaded from https://onlinelibrary.wiley.com/doi/10.1002/prot.22987 by Koc University, Wiley Online Library on [16/02/2026]. See the Terms and Conditions (https://onlinelibrary.wiley.com/terms-and-conditions) on Wiley Online Library for rules of use; OA articles are governed by the applicable Creative Commons License
M.GaoandJ.Skolnick
| TableII |     |     |     |     |     |     | Availability |     |     |     |     |     |     |     |
| ------- | --- | --- | --- | --- | --- | --- | ------------ | --- | --- | --- | --- | --- | --- | --- |
ParametersforCalculatingtheLocationandScaleParametersinEq.(6)
|     |     |     |     |     |     |     | The | data | sets and IS-score |     | software | package | including |     |
| --- | --- | --- | --- | --- | --- | --- | --- | ---- | ----------------- | --- | -------- | ------- | --------- | --- |
IS-score
|     |     |     |     |     |     |     | the source | code | are freely | available | at  | http://cssb.biology. |     |     |
| --- | --- | --- | --- | --- | --- | --- | ---------- | ---- | ---------- | --------- | --- | -------------------- | --- | --- |
gatech.edu/isscore.
| Parameters |     |     | L       | <55 |     | L (cid:6)55 |         |     |     |     |     |     |     |     |
| ---------- | --- | --- | ------- | --- | --- | ----------- | ------- | --- | --- | --- | --- | --- | --- | --- |
|            |     |     | Q       |     |     | Q           |         |     |     |     |     |     |     |     |
| a          |     |     | 0.0806  |     |     | 0.0794      |         |     |     |     |     |     |     |     |
| b          |     |     | 0.0034  |     |     | 0.0033      |         |     |     |     |     |     |     |     |
| c          |     |     | 0.0277  |     |     | 0.0794      | RESULTS |     |     |     |     |     |     |     |
| d          |     |     | 20.0040 |     |     | 20.0054     |         |     |     |     |     |     |     |     |
IS-scoreversusiTM-score
|              |         |     |        |       |           |           | We        | first compare | the        | performance |         | of the     | IS-score | and |
| ------------ | ------- | --- | ------ | ----- | --------- | --------- | --------- | ------------- | ---------- | ----------- | ------- | ---------- | -------- | --- |
|              |         |     |        |       |           |           | iTM-score | on            | evaluating | the         | quality | of docking | models.  |     |
| The criteria | adopted |     | by the | CAPRI | assessors | for model |           |               |            |             |         |            |          |     |
following29: For this comparison, we selected 425 near native and
| quality evaluation |     | are | the |     |     |     |       |           |         |        |      |             |     |       |
| ------------------ | --- | --- | --- | --- | --- | --- | ----- | --------- | ------- | ------ | ---- | ----------- | --- | ----- |
|                    |     |     |     |     |     |     | 5,232 | incorrect | docking | models | from | a Dockgound |     | decoy |
(cid:5) High: f (cid:6) 0.5 & (lRMSD (cid:4) 1 A˚ || iRMSD (cid:4) 1 A˚), set generated with unbound protein structures (see Meth-
nat
(cid:6) > A˚ > ods). As shown in Figure 3(A), the distributions of the
| (cid:5) Medium: | (f         | 0.5   | & lRMSD |       | 1 & iRMSD |              | 1         |     |             |     |           |         |     |        |
| --------------- | ---------- | ----- | ------- | ----- | --------- | ------------ | --------- | --- | ----------- | --- | --------- | ------- | --- | ------ |
|                 | nat        |       |         |       |           |              | IS-scores | for | near native | and | incorrect | docking |     | models |
| A˚) ||          | (f (cid:6) | 0.3 & | f <     | 0.5 & | lRMSD     | (cid:4) 5 A˚ | &         |     |             |     |           |         |     |        |
|                 | nat        |       | nat     |       |           |              |           |     |             |     |           |         |     |        |
iRMSD (cid:4) 2 A˚), are well separated. Near native docking models all have
(cid:6) > A˚ > an IS-score above 0.17, and 97% of the IS-scores >0.25,
| (cid:5) Acceptable: | (f  |     | 0.3 & lRMSD |     | 5 & | iRMSD |     |     |     |     |     |     |        |     |
| ------------------- | --- | --- | ----------- | --- | --- | ----- | --- | --- | --- | --- | --- | --- | ------ | --- |
|                     |     | nat |             |     |     |       |     |     |     |     |     |     | <0.12. |     |
2 A˚) || (f (cid:6) 0.1 & f < 0.3 & lRMSD (cid:4) 10 A˚ & whereas incorrect models all have the scores By
|     | nat |     | nat |     |     |     |     |     |     |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
(cid:4) A˚), comparison, an overlapping regime in the iTM-scores is
| iRMSD              | 4   |       |           |     |         |       |          |         |      |        |     |           |         |     |
| ------------------ | --- | ----- | --------- | --- | ------- | ----- | -------- | ------- | ---- | ------ | --- | --------- | ------- | --- |
|                    |     |       |           |     |         |       | observed | between | near | native | and | incorrect | models. |     |
| (cid:5) Incorrect: | f   | < 0.1 | || (lRMSD | >   | 10 A˚ & | iRMSD | >        |         |      |        |     |           |         |     |
nat
4 A˚). Incorrect docking models have their iTM-scores ranging
|              |               |       |           |         |             |          | from    | 0.33 to    | 0.68; and              | a peak  | is observed |            | at 0.5. | The    |
| ------------ | ------------- | ----- | --------- | ------- | ----------- | -------- | ------- | ---------- | ---------------------- | ------- | ----------- | ---------- | ------- | ------ |
| The notions  |               | & and | || denote | logical | conjunction | and      |         |            |                        |         |             |            |         |        |
|              |               |       |           |         |             |          | peak    | is due     | to the superimposition |         |             | of one     | side    | of the |
| disjunction, | respectively. |       | It        | should  | be noted    | that the |         |            |                        |         |             |            |         |        |
|              |               |       |           |         |             |          | protein | interface. | Most                   | unbound | protein     | structures |         | used   |
A˚
CAPRI assessors employed distance cutoffs of 5 and 10 for docking are structurally very close to their bound
| to define | interfacial | residues |     | separately | for calculating | f   |            |        |          |        |          |      |      |        |
| --------- | ----------- | -------- | --- | ---------- | --------------- | --- | ---------- | ------ | -------- | ------ | -------- | ---- | ---- | ------ |
|           |             |          |     |            |                 | nat | structural | forms. | In these | cases, | at least | half | of a | native |
and iRMSD. In this study, we only used the final assess- interface can be superimposed to its counterpart in
ments (i.e., High, Medium, Acceptable, and Incorrect) a docking model, despite the fact that the other side of
| provided | by the | CAPRI | assessors. |     |     |     |               |     |             |      |           |          |     |       |
| -------- | ------ | ----- | ---------- | --- | --- | --- | ------------- | --- | ----------- | ---- | --------- | -------- | --- | ----- |
|          |        |       |            |     |     |     | the interface |     | is far away | from | it native | position |     | in an |
Figure 3
ComparisonoftheiTM/IS-scoresforassessingthequalityofproteindockingmodels.(A)Scoredistributionsofincorrectdockingmodelsandof
nearnativedockingmodels.iTMSandISSdenoteiTM-scoreandISSscore,respectively.(B)ROCcurvesofsensitivityversusfalse-positiverate.
[Colorfigurecanbeviewedintheonlineissue,whichisavailableatwileyonlinelibrary.com.]
1628 PROTEINS

AssessmentofProteinDockingModels
Figure 4
Qualityassessmentsof152,600dockingmodelsgeneratedfor1,526proteincomplexes.(A)NumberofdockingmodelsaccordingtotheIS-scoreP
values.Boxplotsofdockingmodelsaccordingto(B)fractionofnativecontactspreservedinmodels,(C)interfacial,and(D)ligandRMSDs.The
lower,middleandupperquartilesofeachboxarethe25th,50th,and75thpercentile;whiskersextendtoadistanceofupto1.98timesthe
interquartilerange.Outliersandmeansarerepresentedbycircles.[Colorfigurecanbeviewedintheonlineissue,whichisavailableat
wileyonlinelibrary.com.]
incorrect model. Such superimposition gives a significant near-native models, and the false-positive rate is the frac-
iTM-score >0.4, as overlapping the score regime of the tion of incorrect models. The ROC curves were obtained
near native models from 0.4 to 0.9. by varying the thresholds of the iTM/IS-score. The IS-
The performance of IS-score and iTM-score is further score has a perfect ROC curve with the value of AUC 0.2
displayed in the Receiver Operating Characteristic (ROC) (Area Under Curve up to a 20% false-positive rate) of 1,
curves [Fig. 3(B)], where the sensitivity is the fraction of whereas the iTM-score has an AUC value of 0.76.
0.2
PROTEINS 1629
10970134,
2011,
5,
Downloaded
from
https://onlinelibrary.wiley.com/doi/10.1002/prot.22987
by
Koc
University,
Wiley
Online
Library
on
[16/02/2026].
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

 10970134, 2011, 5, Downloaded from https://onlinelibrary.wiley.com/doi/10.1002/prot.22987 by Koc University, Wiley Online Library on [16/02/2026]. See the Terms and Conditions (https://onlinelibrary.wiley.com/terms-and-conditions) on Wiley Online Library for rules of use; OA articles are governed by the applicable Creative Commons License
M.GaoandJ.Skolnick
Figure 5
Twodockingmodelsfor(A)aputativecitratelyase(PDBcode:1xr4,chainAandB)and(B)anaminotransferase(PDBcode:1dty,chainAand
B).Ineachsnapshot,thetwochainsfromdockingmodelarecoloredincyan/orange,andthecorrespondingchainsinthenativestructuresare
coloredinblue/red.Forclarity,interface/noninterfaceregionsareshowninsolid/transparentcolors,respectively.Overlappedinterfaceregionsare
indicatedbyagreenbackgroundin(B).MolecularimageswerecreatedwithVMD.39Thelengthoftheinterfaceandthenumberofinterfacial
| contactsaredenotedasN |     | andN | ,respectively. |     |     |     |     |     |     |     |     |     |     |
| --------------------- | --- | ---- | -------------- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
|                       |     | res  | con            |     |     |     |     |     |     |     |     |     |     |
Overall, the analysis demonstrates that a similarity metric the native structure, validating the estimated high P value
|     |     |     |     |     |     |     | 3   | 223. |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | ---- | --- | --- | --- | --- | --- |
based purely on geometric distances has an intrinsic flaw of 6.5 10 In Figure 5(B), the docking model has a
27,
for evaluating docking models and that the IS-score P of 5.4 3 10 due to the maintenance of 22 native
yields a much more accurate assessment by taking inter- contacts, despite a different orientation from the native
| facial contacts             | explicitly |     | into account. |     |     |     | docking   | pose. |                           |        |        |      |        |
| --------------------------- | ---------- | --- | ------------- | --- | --- | --- | --------- | ----- | ------------------------- | ------ | ------ | ---- | ------ |
|                             |            |     |               |     |     |     | Virtually | all   | insignificant             | models | at P > | 0.01 | has an |
|                             |            |     |               |     |     |     | iRMSD>3A˚ |       | <10%.About1%ofdockingmod- |        |        |      |        |
| Discriminatingdockingmodels |            |     |               |     |     |     |           | andf  |                           |        |        |      |        |
nat
|            |             |                |               |              |            |              | elsexhibitaninterfacethatbearsasignificant |            |           |             |               | similarityto |     |
| ---------- | ----------- | -------------- | ------------- | ------------ | ---------- | ------------ | ------------------------------------------ | ---------- | --------- | ----------- | ------------- | ------------ | --- |
| To further | examine     | whether        |               | the IS-score | returns    | a rea-       |                                            |            |           |             |               | 3            | 26. |
|            |             |                |               |              |            |              | the native                                 | interface  | with      | a P between | 0.01          | and 1        | 10  |
| sonable    | estimate    | of statistical | significance, |              | we         | further per- |                                            |            |           |             |               |              |     |
|            |             |                |               |              |            |              | These model                                | interfaces | typically |             | have a iRMSD  | between      | 5   |
| formed     | large scale | tests          | on a          | total        | of 152,600 | docking      |                                            |            |           |             |               |              |     |
|            |             |                |               |              |            |              | and 10 A˚                                  | and        | preserve  | 10% to      | 30% of native | contacts.    |     |
| models     | for the     | 1,526 target   | complexes.    |              | Each       | model was    |                                            |            |           |             |               |              |     |
Theyusuallyoverlapapartofthenativeinterface.
| assessed          | according | to the       | IS-score      | with     | respect  | to the   |     |     |     |     |     |     |     |
| ----------------- | --------- | ------------ | ------------- | -------- | -------- | -------- | --- | --- | --- | --- | --- | --- | --- |
| native structure. |           | As expected, |               | the vast | majority | (96%) of |     |     |     |     |     |     |     |
| these models      | have      | an           | insignificant |          | IS-score | with P   | >   |     |     |     |     |     |     |
AssessingCAPRImodels
| 0.01, while | a small | fraction | (3.2%) |     | of docking | models |     |     |     |     |     |     |     |
| ----------- | ------- | -------- | ------ | --- | ---------- | ------ | --- | --- | --- | --- | --- | --- | --- |
resemble the native structure at a high level of similarity Finally, we applied the IS-score to assess the quality of
210
with P < 1 3 10 [Fig. 4(A)]. docking models submitted by various research groups for
As shown in Figure 4(B,C), all docking models within 10 recent CAPRI targets. The results of the IS-score eval-
A˚
2.5 iRMSD from native structures or with a f value uations are compared to the official assessments provided
nat
>30% have a significant P better than 1 3 10 26, mostly, by the CAPRI organizers, who categorized each model
3 210.
better than 1 10 Conversely, almost all interfaces into one of four groups: Incorrect, Acceptable, Medium,
26
with P < 1 3 10 have an iRMSD of less than 2.5 A˚ and High, according to iRMSD, lRMSD, and f (see
nat
and a f nat of more than 30%. In rare exceptions, a dock- Methods). A total of 2,874 Incorrect, 117 Acceptable, 59
|     |     |     | <   | 3   | 26, |     |     |     |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
ing model has a significant P 1 10 while exhibit- Medium, and 16 High quality models for these ten tar-
ing a relatively high iRMSD/lRMSD >3/8 A˚ and low gets were evaluated. Consistent with the CAPRI assess-
<30%.
native contacts These cases are from docking very ments, the overall distributions of the four groups of
large complexes with usually more than 150 interfacial docking models are clearly separated according to either
amino acids. Two cases are shown in Figure 5. Despite a the IS-scores or their P values (Fig. 6). The means of the
|     |     | A˚, |     |     |     |     |     |     | 0.08/20.21 |     | 0.26/25.7 |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | ---------- | --- | --------- | --- | --- |
high lRMSD of 9.5 visual inspection suggests that the IS-scores/Log P are (I), (A), 0.48/
10
docking model shown in Figure 5(A) resembles very well 214.0 (M), and 0.69/221.2 (H), respectively.
1630 PROTEINS

AssessmentofProteinDockingModels
Figure 6
DistributionofCAPRImodelsaccordingto(A)theIS-scorePvaluesand(B)theIS-score.LegendsindicatemodelqualityprovidedbytheCAPRI
assessors.(C)Oneexample(modelID:T26_P41.M02)ofCAPRIdockingmodelsfortargetT26.ThemodelwascategorizedasIncorrectaccording
totheCAPRIassessors,butshowssignificantinterfacesimilarity.ThecoloringschemeisthesameasthatemployedinFigure5.
Out of 192 models with better or acceptable quality, attributed to two main reasons. First, the CAPRI assess-
174 (91%) and 185 (96%) have a significant P < 0.01 ment uses f and RMSDs, with a size-dependence issue,
nat
and 0.05, respectively. Only seven Acceptable models whereas the IS-score takes the length effect into account.
have a P > 0.05. These models, from targets T24, T25, Second, the IS-score only considers interface similarity
T27, and T29, have about 10% to 15% native contacts but ignores global orientation. A slight rotation could
correctly modeled. However, the numbers of preserved lead to a large lRMSD, despite the fact that the iRMSD is
native contacts are small, considering that their native relatively small. An example of an Incorrect model with
interfaces consist of about 40 native contacts or less. On significant interface similarity is shown in Figure 6(C).
the other hand, a total of 263 models have a significant Visual inspection suggests that the model has good inter-
similarity to their target interface at a P < 0.01 or better. face similarity at iRMSD of 4.3 A˚ and 27% f . However, nat
Among these significant models, 89 were assigned as a slight tilt around the interface leads to a large lRMSD
Incorrect. Most of these significant Incorrect models are of 12 A˚.
from targets T26 and T32, which have relatively large Figure 7 shows the quality of individual docking mod-
interfaces with more than 60 native contacts. The differ- els for each target. For all targets with the exception of
ence between the CAPRI and IS-score assessments can be T25, unbound or homology structures were provided as
PROTEINS 1631
10970134,
2011,
5,
Downloaded
from
https://onlinelibrary.wiley.com/doi/10.1002/prot.22987
by
Koc
University,
Wiley
Online
Library
on
[16/02/2026].
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

 10970134, 2011, 5, Downloaded from https://onlinelibrary.wiley.com/doi/10.1002/prot.22987 by Koc University, Wiley Online Library on [16/02/2026]. See the Terms and Conditions (https://onlinelibrary.wiley.com/terms-and-conditions) on Wiley Online Library for rules of use; OA articles are governed by the applicable Creative Commons License
M.GaoandJ.Skolnick
Figure 7
Thequalityofindividualdockingmodelssubmittedbydifferentresearchgroupsfor10CAPRItargets.ThetargetIDisshowninthelowerleft
cornerofeachplot.ThehorizontaldashedlinesarelocatedatP50.05accordingtotheIS-score.[Colorfigurecanbeviewedintheonlineissue,
whichisavailableatwileyonlinelibrary.com.]
the starting structures for docking experiments. Overall, was found for T28, the top ranked docking model of the
it is clear that higher ranked models have better quality, same target has a marginal P value of 0.047.
| consistent  | with the official    | assessments. In             | particular, for |                  |           |             |          |
| ----------- | -------------------- | --------------------------- | --------------- | ---------------- | --------- | ----------- | -------- |
| targets T35 | and T36, where       | only one Acceptable         | model           |                  |           |             |          |
|             |                      |                             |                 | DISCUSSION       | AND       | CONCLUSION  |          |
| was found,  | these two Acceptable | models were                 | the best as     |                  |           |             |          |
| assessed    | by the P-value of    | the IS-score. Additionally, | cor-            |                  |           |             |          |
|             |                      |                             |                 | Currently, iRMSD | and lRMSD | are metrics | commonly |
responding to the assessment that no Acceptable model employed in docking studies. The major advantages of
Figure 8
AnexampleillustratesthatlocalstructuralsimilarityiscapturedbytheIS-scorebutnotbyiRMSD.Themodelstructureissuperimposedontothe
nativestructurein(A)topviewand(B)sideview.Themodelstructureoverlapsthenativestructure(PDBcode2cwq)intheinterfaceregion,
exceptfortwohelicalsegments(labeledasH1andH2)exhibitinganalmost180degreerotation,oneofwhichisindicatedbyagreyarrow.
Interfacialregionsofthenativestructureandthecorrespondingresiduesinthemodelstructureareshowninsolidcolors,andotherregionsare
transparent.
1632 PROTEINS

 10970134, 2011, 5, Downloaded from https://onlinelibrary.wiley.com/doi/10.1002/prot.22987 by Koc University, Wiley Online Library on [16/02/2026]. See the Terms and Conditions (https://onlinelibrary.wiley.com/terms-and-conditions) on Wiley Online Library for rules of use; OA articles are governed by the applicable Creative Commons License
AssessmentofProteinDockingModels
the RMSD metrics are twofold: first, the overall quality randomly selected docking models. Virtually all highly
26
of a docking model is guaranteed if one uses a very con- significant interfaces with P > 10 are native-like, and
servative RMSD criterion; second, the calculation of conversely, all native-like docking models display a highly
|     |     |     |     |     |     |     |     |     | >   | 26, |     | >10 210. |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | -------- | --- | --- |
RMSD is very straightforward. However, RMSD metrics significant P 10 mostly By contrast, insig-
also have two significant disadvantages. First, it is well nificant models with P > 10 26 have a iRMSD >3 A˚ and
|     |     |     |     |     |     |     |     | <   |     |     |     |     |     | 3 26 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | ---- |
known that the statistical significance of a given RMSD a f 10%. Models with P between 0.01 and 1 10
nat
value is length dependent [e.g., Fig. 1(A)]. As a result, have some interfacial similarity, but may exhibit a rota-
there is no straightforward relationship between RMSD tion that gives a relatively large lRMSD.
values and the statistical significance of docking models. The IS-score is further applied to evaluate the docking
This is reflected in a simple fact that, at the same iRMSD models for ten recent CAPRI targets. Overall, the evalua-
A˚),
value (e.g., 3 to build a docking model for a 100-resi- tion of the IS-score is consistent with the official CAPRI
due interface is more difficult than for a 30-residue inter- assessment. On average, the mean of the IS-scores are
face. In addition, RMSD metrics are global metrics, 0.26, 0.48, and 0.69, for Acceptable, Medium, and High
meaning that local similarity may not be properly charac- resolution models, respectively. However, it appears that
terized by RMSD. One extreme example is shown in Fig- the official assessment is somewhat conservative. Accord-
ure 8, where the docking model has a highly significant ing to the P values of the IS-scores, we identified quite a
218),
IS-score of 0.43 (P 5 2 3 10 despite a large iRMSD few models whose significance is underestimated. The IS-
value of 13.2 A˚, caused by the rotations of two helical score scheme is conceptually simple and statistically
segments. Other than the two helical segments, the re- sound. One further application of the scheme is to use it
mainder (60%) of the interface superimposes with an as the objective function for method optimization.
A˚
| RMSD of    | less than  | 2   | between | the model |         | and the | native |     |     |     |     |     |     |     |
| ---------- | ---------- | --- | ------- | --------- | ------- | ------- | ------ | --- | --- | --- | --- | --- | --- | --- |
| structure. | Obviously, | the | model   | in this   | example | is      | not a  |     |     |     |     |     |     |     |
REFERENCES
| random        | prediction. | For          | the purpose |                  | of assessing |            | a dock- |              |         |         |                    |     |             |     |
| ------------- | ----------- | ------------ | ----------- | ---------------- | ------------ | ---------- | ------- | ------------ | ------- | ------- | ------------------ | --- | ----------- | --- |
| ing method,   | it          | is important |             | to differentiate |              | such       | a case  |              |         |         |                    |     |             |     |
|               |             |              |             |                  |              |            |         | 1. Aebersold | R, Mann | M. Mass | spectrometry-based |     | proteomics. | Na- |
| from a random |             | model        | prediction. | Overall,         |              | one should | be      |              |         |         |                    |     |             |     |
ture2003;422:198–207.
cautious in using RMSD metrics to assess the quality of 2. RigautG,ShevchenkoA,RutzB, WilmM,MannM,SeraphinB.A
a docking model. generic protein purification method for protein complex characteri-
We have introduced and examined the performance of zationandproteomeexploration.NatBiotechnol1999;17:1030–1032.
|             |             |             |           |                    |             |           |     | 3. Uetz P,                             | Giot L, Cagney |           | G, Mansfield   | TA,   | Judson RS,        | Knight JR, |
| ----------- | ----------- | ----------- | --------- | ------------------ | ----------- | --------- | --- | -------------------------------------- | -------------- | --------- | -------------- | ----- | ----------------- | ---------- |
| two scoring | schemes,    | the         | iTM-score |                    | and the     | IS-score, | for |                                        |                |           |                |       |                   |            |
|             |             |             |           |                    |             |           |     | LockshonD,NarayanV,SrinivasanM,Pochart |                |           |                |       | P,Qureshi-EmiliA, |            |
| use in      | assessing   | the quality |           | of protein-protein |             | docking   |     |                                        |                |           |                |       |                   |            |
|             |             |             |           |                    |             |           |     | Li Y, Godwin                           | B,             | Conover   | D, Kalbfleisch | T,    | Vijayadamodar     | G, Yang    |
| models.     | Both scores | are         | able to   | detect             | significant | substruc- |     |                                        |                |           |                |       |                   |            |
|             |             |             |           |                    |             |           |     | MJ, Johnston                           | M,             | Fields S, | Rothberg       | JM. A | comprehensive     | analysis   |
ture similarity if it exists. While the iTM-score is based of protein-protein interactions in Saccharomyces cerevisiae. Nature
| on geometric    | distances, |               | the IS-score |            | combines   | both         | inter- | 2000;403:623–627. |            |           |               |            |          |             |
| --------------- | ---------- | ------------- | ------------ | ---------- | ---------- | ------------ | ------ | ----------------- | ---------- | --------- | ------------- | ---------- | -------- | ----------- |
|                 |            |               |              |            |            |              |        | 4. Gavin AC,      | Aloy       | P, Grandi | P, Krause     | R, Boesche | M,       | Marzioch M, |
| facial contacts |            | and geometric |              | distances. |            | In benchmark |        |                   |            |           |               |            |          |             |
|                 |            |               |              |            |            |              |        | Rau C, Jensen     | LJ,        | Bastuck   | S, Dumpelfeld | B,         | Edelmann | A, Heurtier |
| tests of        | 425 near   | native        | models       | and        | 5,232      | randomly     |        |                   |            |           |               |            |          |             |
|                 |            |               |              |            |            |              |        | MA, Hoffman       | V,         | Hoefert   | C, Klein      | K, Hudak   | M,       | Michon AM,  |
| related,        | incorrect  | models,       | generated    | from       | rigid-body |              | dock-  |                   |            |           |               |            |          |             |
|                 |            |               |              |            |            |              |        | Schelder          | M, Schirle | M,        | Remor         | M, Rudi T, | Hooper   | S, Bauer A, |
ing of unbound protein structures, the IS-score achieves a BouwmeesterT, CasariG,DrewesG,NeubauerG,RickJM,Kuster
perfect classification at an AUC value of 1, whereas the B, Bork P, Russell RB, Superti-Furga G. Proteome survey reveals
0.2
modularityoftheyeastcellmachinery.Nature2006;440:631–636.
| iTM-score | gives | an inferior |     | performance |     | at an | AUC 0.2 |     |     |     |     |     |     |     |
| --------- | ----- | ----------- | --- | ----------- | --- | ----- | ------- | --- | --- | --- | --- | --- | --- | --- |
5. GiotL,BaderJS,BrouwerC,ChaudhuriA,KuangB,LiY,HaoYL,
| value of        | 0.76. The | main    | issue          | with the | iTM-score |               | is that |         |        |           |                  |     |            |           |
| --------------- | --------- | ------- | -------------- | -------- | --------- | ------------- | ------- | ------- | ------ | --------- | ---------------- | --- | ---------- | --------- |
|                 |           |         |                |          |           |               |         | Ooi CE, | Godwin | B, Vitols | E, Vijayadamodar |     | G, Pochart | P, Machi- |
| the interaction |           | pose is | not explicitly |          | taken     | into account. |         |         |        |           |                  |     |            |           |
neniH,WelshM,KongY,ZerhusenB,MalcolmR,VarroneZ,Col-
As a result, an artificially high iTM-score may be obtained lisA,MintoM,BurgessS,McDanielL,StimpsonE,SpriggsF,Wil-
through the superimposition of one side of the interface, liamsJ,NeurathK,IoimeN,AgeeM,VossE,FurtakK,RenzulliR,
|           |       |         |               |     |     |          |      | Aanensen | N, Carrolla | S,  | Bickelhaupt | E, Lazovatsky |     | Y, DaSilva A, |
| --------- | ----- | ------- | ------------- | --- | --- | -------- | ---- | -------- | ----------- | --- | ----------- | ------------- | --- | ------------- |
| while the | other | side of | the interface | may | be  | far away | from |          |             |     |             |               |     |               |
ZhongJ,StanyonCA,FinleyRL,WhiteKP,BravermanM,JarvieT,
| its native | position. | The       | issue     | is intrinsic |     | to all     | scoring |         |          |        |             |     |         |              |
| ---------- | --------- | --------- | --------- | ------------ | --- | ---------- | ------- | ------- | -------- | ------ | ----------- | --- | ------- | ------------ |
|            |           |           |           |              |     |            |         | Gold S, | Leach M, | Knight | J, Shimkets | RA, | McKenna | MP, Chant J, |
| functions  | based     | solely on | geometric | distances.   |     | By compar- |         |         |          |        |             |     |         |              |
RothbergJM.AproteininteractionmapofDrosophilamelanogaster.
ison, the introduction of the contact overlap factor in the Science2003;302:1727–1736.
IS-score scheme eliminates this issue. Since the IS-score is 6. KroganNJ,CagneyG,YuHY,ZhongGQ,GuoXH,IgnatchenkoA,
|           |         |       |           |     |          |             |     | Li J, Pu | SY, Datta | N, Tikuisis | AP, | Punna T, | Peregrin-Alvarez | JM, |
| --------- | ------- | ----- | --------- | --- | -------- | ----------- | --- | -------- | --------- | ----------- | --- | -------- | ---------------- | --- |
| dependent | on side | chain | contacts, | it  | requires | an accurate |     |          |           |             |     |          |                  |     |
ShalesM,ZhangX,DaveyM,RobinsonMD,PaccanaroA,BrayJE,
| side-chain  | reconstruction |                | procedure |           | in order | to evaluate |     |        |            |             |     |          |          |            |
| ----------- | -------------- | -------------- | --------- | --------- | -------- | ----------- | --- | ------ | ---------- | ----------- | --- | -------- | -------- | ---------- |
|             |                |                |           |           |          |             |     | Sheung | A, Beattie | B, Richards | DP, | Canadien | V, Lalev | A, Mena F, |
| the quality | of a           | coarse-grained |           | Ca model. |          |             |     |        |            |             |     |          |          |            |
WongP,StarostineA,CaneteMM,VlasblomJ,WuS,OrsiC,Col-
For a proper model quality assessment, it is important lins SR, Chandran S, Haw R, Rilstone JJ, Gandi K, Thompson NJ,
to assess the statistical significance of predicted models. Musso G, St Onge P, Ghanny S, Lam MHY, Butland G, Altaf-Ui
Using random interfaces as the background, we have AM, Kanaya S, Shilatifard A, O’Shea E, Weissman JS, Ingles CJ,
HughesTR,ParkinsonJ,GersteinM,WodakSJ,EmiliA,Greenblatt
| derived          | statistical | models         | for | estimating   | the | significance |         |            |           |            |           |     |           |             |
| ---------------- | ----------- | -------------- | --- | ------------ | --- | ------------ | ------- | ---------- | --------- | ---------- | --------- | --- | --------- | ----------- |
|                  |             |                |     |              |     |              |         | JF. Global | landscape | of protein | complexes | in  | the yeast | Saccharomy- |
| of the IS-score. |             | The estimation |     | is validated |     | on           | 156,200 |            |           |            |           |     |           |             |
cescerevisiae.Nature2006;440:637–643.
|     |     |     |     |     |     |     |     |     |     |     |     |     | PROTEINS | 1633 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | -------- | ---- |

M.GaoandJ.Skolnick
7. Li SM, Armstrong CM, Bertin N, Ge H, Milstein S, Boxem M, 22. Wang C, Bradley P, Baker D. Protein-protein docking with back-
Vidalain PO, Han JDJ, Chesneau A, Hao T, Goldberg DS, Li N, boneflexibility.JMolBiol2007;373:503–519.
MartinezM,RualJF,LameschP,XuL,TewariM,WongSL,Zhang 23. Gabb HA, Jackson RM, Sternberg MJE. Modelling protein docking
LV, Berriz GF, Jacotot L, Vaglio P, Reboul J, Hirozane-Kishikawa T, using shape complementarity, electrostatics and biochemical infor-
LiQR,GabelHW,ElewaA,BaumgartnerB,RoseDJ,YuHY,Bosak mation.JMolBiol1997;272:106–120.
S,SequerraR,FraserA,MangoSE,SaxtonWM,StromeS,vanden 24. Gao M, Skolnick J. Structural space of protein-protein interfaces is
Heuvel S, Piano F, Vandenhaute J, Sardet C, Gerstein M, Doucet- degenerate, close to complete, and highly connected. Proc Natl
te-Stamm L, Gunsalus KC, Harper JW, Cusick ME, Roth FP, Hill AcadSciUSA2010;107:22517–22522.
DE, Vidal M. A map of the interactome networkof the metazoan 25. Bonvin AM. Flexible protein-protein docking. Curr Opin Struct
C.elegans.Science2004;303:540–543. Biol2006;16:194–200.
8. Russell RB, Alber F, Aloy P, Davis FP, Korkin D, Pichaud M, Topf 26. Kastritis PL, Bonvin AM. Are scoring functions in protein-protein
M, Sali A. A structural perspective on protein-protein interactions. docking ready to predict interactomes? Clues from a novel binding
CurrOpinStructBiol2004;14:313–324. affinitybenchmark.JProteomeRes2010;9:2216–2225.
9. Tuncbag N, Gursoy A, Guney E, Nussinov R, Keskin O. Architec- 27. Janin J. Assessing predictions of protein-protein interaction: the
tures and functional coverage of protein-protein interfaces. J Mol CAPRIexperiment.ProteinSci2005;14:278–283.
Biol2008;381:785–802. 28. Lensink MF, Mendez R, Wodak SJ. Docking and scoring protein
10. AloyP,RussellRB.Interrogatingproteininteractionnetworksthrough complexes:CAPRI3rdedition.Proteins2007;69:704–718.
structuralbiology.ProcNatlAcadSciUSA2002;99:5896–5901. 29. Mendez R, Leplae R, De Maria L, Wodak SJ. Assessment of blind
11. Chen HL, Skolnick J. M-TASSER: an algorithm for protein quater- predictions of protein-protein interactions: current status of dock-
narystructureprediction.BiophysJ2008;94:918–928. ingmethods.Proteins2003;52:51–67.
12. LuL,LuH,SkolnickJ.MULTIPROSPECTOR:analgorithmforthe 30. Reva BA, Finkelstein AV, Skolnick J. What is the probability of a
prediction of protein-protein interactions by multimeric threading. chance prediction of a protein structure with an rmsd of 6 ang-
Proteins2002;49:350–364. strom?FoldDes1998;3:141–147.
13. Gunther S, May P, Hoppe A, Frommel C, Preissner R. Docking 31. Zemla A. LGA: a method for finding 3D similarities in protein
withoutdocking:ISEARCH–predictionofinteractionsusingknown structures.NucleicAcidsRes2003;31:3370–3374.
interfaces.Proteins2007;69:839–844. 32. Siew N, Elofsson A, Rychiewski L, Fischer D. MaxSub: an auto-
14. Keskin O, Nussinov R, Gursoy A. PRISM: protein-protein interaction mated measure for the assessment of protein structure prediction
predictionbystructuralmatching.MethodsMolBiol2008;484:505–521. quality.Bioinformatics2000;16:776–785.
15. Sinha R, Kundrotas PJ, Vakser IA. Docking by structural similarity 33. Zhang Y, Skolnick J. Scoring function for automated assessment of
atprotein-proteininterfaces.Proteins2010;78:3235–3241. proteinstructuretemplatequality.Proteins2004;57:702–710.
16. Sandak B, Wolfson HJ, Nussinov R. Flexible docking allowing 34. Levitt M, Gerstein M. A unified statistical framework for sequence
induced fit in proteins: insights from an open to closed conforma- comparison and structure comparison. Proc Natl Acad Sci USA
tionalisomers.Proteins1998;32:159–174. 1998;95:5913–5920.
17. Chen R, Li L, Weng ZP. ZDOCK: an initial-stage protein-docking 35. OrtizAR,StraussCEM,OlmeaO.MAMMOTH(Matchingmolecu-
algorithm.Proteins2003;52:80–87. larmodelsobtainedfromtheory):anautomatedmethodformodel
18. DominguezC,BoelensR,BonvinA.HADDOCK:aprotein-protein comparison.ProteinSci2002;11:2606–2621.
dockingapproachbasedonbiochemicalorbiophysicalinformation. 36. Gao M, Skolnick J. iAlign: a method for the structural compari-
JAmChemSoc2003;125:1731–1737. son of protein-protein interfaces. Bioinformatics 2010;26:2259–
19. Fernandez-Recio J, Totrov M, Abagyan R. Soft protein-protein 2265.
dockingininternalcoordinates.ProteinSci2002;11:280–291. 37. Kabsch W. Solution for best rotation to relate two sets of vectors.
20. Kozakov D, Brenke R, Comeau SR, Vajda S. PIPER: an FFT-based ActaCrystallogrA1976;32:922–923.
protein docking program with pairwise potentials. Proteins 38. Liu S, Gao Y, Vakser IA. DOCKGROUND protein-protein docking
2006;65:392–406. decoyset.Bioinformatics2008;24:2634–2635.
21. Vakser IA. Protein docking for low-resolution structures. Protein 39. HumphreyW,DalkeA,SchultenK.VMD:visualmoleculardynam-
Eng1995;8:371–377. ics.JMolGraphics1996;14:33–38.
1634 PROTEINS
10970134,
2011,
5,
Downloaded
from
https://onlinelibrary.wiley.com/doi/10.1002/prot.22987
by
Koc
University,
Wiley
Online
Library
on
[16/02/2026].
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