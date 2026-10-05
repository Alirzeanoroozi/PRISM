Bioinformatics,38(2),2022,554–555
doi:10.1093/bioinformatics/btab596
AdvanceAccessPublicationDate:9September2021
ApplicationsNote
Structural bioinformatics
INTERCAAT: identifying interface residues between
macromolecules
Steven Grudman, J. Eduardo Fajardo andAndras Fiser *
DepartmentofSystemsandComputationalBiology,AlbertEinsteinCollegeofMedicine,Bronx,NY10461,USA
*Towhomcorrespondenceshouldbeaddressed.
AssociateEditor:LenoreCowen
ReceivedonJune9,2021;editorialdecisiononAugust12,2021;revisedonJuly21,2021
Abstract
Summary: The Interface Contact definition with Adaptable Atom Types (INTERCAAT) was developed to determine
theatomicinteractionsbetweenmoleculesthatformaknownthreedimensionalstructure.First,INTERCAATcreates
a Voronoi tessellation where each atom acts as a seed. Interactions are defined by atoms that share a hyperplane
andwhosedistanceislessthanthesumofeachatoms’VanderWaalsradiiplusthediameterofasolventmolecule.
Interactingatomsarethenclassifiedandinteractionsarefilteredbasedoncompatibility.INTERCAATimplementsan
adaptiveatomclassificationmethod;therefore,itcanexploreinterfacesbetweenavarietymacromolecules.
Availabilityandimplementation:Sourcecodeisfreelyavailableat:https://gitlab.com/fiserlab.org/intercaat.
Contact:andras.fiser@einsteinmed.org
Supplementaryinformation:SupplementarydataareavailableatBioinformaticsonline.
1Introduction considering the Van der Waals radii among heavy atoms are not
drastically different. Atomic interactions can be further filtered to
Exploring interfaces of macromolecular interactions from Protein show only ‘legitimate’ interactions. Legitimacy depends on the
DataBank(PDB)coordinatefiles(Bermanetal.,2000)isanessen- hydrophobic/hydrophilic properties of the interacting atoms
tial everyday task in bioinformatics. Several software tools have (Sobolevetal.,1999).Atomscanbelongtooneofeightclassesand
beendevelopedtoutilizePDBcoordinatestovisualizeandanalyze ifeachatomclassiscompatible,theirinteractionisconsidered‘legit-
inter and intra molecular interactions (Sobolev et al., 1999). imate’.Foraresidueonthequerychaintobeconsideredaspartof
Determiningresiduesthatformtheinterfacebetweenproteinsissur- theinterface,itis requiredtohavea minimum numberofinterac-
prisinglycomplicated.Arecentstudydemonstratedthatonaverage tions with the interacting chain(s) to prevent accidental classifica-
only about 80% of residue overlap between any two alternative tions. Voronoi tessellations were first used in a protein context in
interfacepredictionmethods(GilandFiser,2019).Thisisduetothe
1974(Richards,1974)buthavesincebeenusedtoinvestigateaser-
subjective definitions guiding these methods, some of which focus
iesofproteinrelatedissuesincludingresiduevolumes,packing,fold-
on changes in solvent accessibility, while others focus on variable
ingandbinding(Poupon,2004).
distance thresholds requiring specific contacts between interacting
residues.Ourcurrenteffortfocusedonestablishingagenericmethod
to accurately assess interfaces using an advanced geometrical ap- 2Softwaredesign
proach,consideringthecompatibilityofinteractions,andproviding
adjustable options that the user can modify to explore alternative INTERCAAT was developed in a Linux environment. It requires
definitions. Another advantage of INTERCAAT is that it uses an threeinputswhilesixadditionalswitchesareoptional.Therequired
adaptive atom classification function and therefore can explore inputs include the name of a PDB file, the chain ID of the query
interactions between a variety of molecules beyond proteins e.g. chainwhoseinterfaceneedstobedetermined,andthechainID(s)of
interactionswithnucleicacidsorlipids. thechain(s)interactingwiththequerychain.Theoptionalswitches
First,INTERCAATparsesaPDBfileandcreatesaVoronoites- includesettingtheminimumrequirednumberofinteractionsofthe
sellation between atoms. A Voronoi diagram is computed via query chain required, whether to display an interaction matrix,
Delaunay triangulation. An atomic interaction is established be- whether to consider class compatibility of interactions, setting the
tweentwoatomsiftheyshareaboundedhyperplaneandarewithin solvent molecule radius, whether to include chains other than the
adistancelessthanthesumoftheatomsVanderWaalsradiiplus queryandinteractingchainsintheVoronoicalculation,andfinally,
the diameter of a solvent molecule. We should point out that our anoptionalfilepathofthePDBfile.Eachoptionalswitchhasde-
methodologytreatstheatomsaspointsintheVoronoitessellation faultvalues;which,alongwithaninputexample,canbedisplayed
andthenasspheresforthedistancecutoff.Thisisasmallconflict withtheprogramshelpfunction.Theoutputdisplayseveryatomic
VCTheAuthor(s)2021.PublishedbyOxfordUniversityPress.Allrightsreserved.Forpermissions,pleasee-mail:journals.permissions@oup.com 554

| INTERCAAT |     |     |     |     |     |                   |     |                |     |                 |           |               | 555      |
| --------- | --- | --- | --- | --- | --- | ----------------- | --- | -------------- | --- | --------------- | --------- | ------------- | -------- |
|           |     |     |     |     |     | these comparisons |     | focused        | on  | protein–protein |           | interactions, | we       |
|           |     |     |     |     |     | should point      | out | that INTERCAAT |     | was             | developed | with          | adaptive |
atomclassificationcapabilitiesanditisnotrestrictedtoprotein–pro-
teininteractions(SupplementaryMaterial).
Thegoalofthecomparisonswastoevaluatetheperformanceof
INTERCAATaswellastodeterminetheoptimalinputforthemin-
|     |     |     |     |     |     | imum atomic | interactions |     | necessary | for a | residue | to be considered |     |
| --- | --- | --- | --- | --- | --- | ----------- | ------------ | --- | --------- | ----- | ------- | ---------------- | --- |
partoftheinterface.Thiswasdonebycomparingthecommoninter-
|     |     |     |     |     |     | face residues | predicted | by       | both INTERCAAT   |     | and         | the other | data-    |
| --- | --- | --- | --- | --- | --- | ------------- | --------- | -------- | ---------------- | --- | ----------- | --------- | -------- |
|     |     |     |     |     |     | bases, the    | unique    | residues | INTERCAAT        |     | predicted   | and the   | unique   |
|     |     |     |     |     |     | residues      | predicted | by the   | other databases. |     | In addition | to        | plotting |
boxandwhiskerplots,wherethewhiskershaveacutoffatthe5th
|     |     |     |     |     |     | and 95th | percentiles, | F   | scores were | calculated | to  | quantify | these |
| --- | --- | --- | --- | --- | --- | -------- | ------------ | --- | ----------- | ---------- | --- | -------- | ----- |
1
|     |     |     |     |     |     | results into | a single | score. | F scores | represent | a tests | accuracy | by  |
| --- | --- | --- | --- | --- | --- | ------------ | -------- | ------ | -------- | --------- | ------- | -------- | --- |
1
measuringtheharmonicmeanutilizingatest’sprecisionandrecall
| Fig.1.ComparisonofINTERCAATagainstthreeotherdatabasesforprotein–pro- |     |     |     |     |     | (Fig.1). |     |     |     |     |     |     |     |
| -------------------------------------------------------------------- | --- | --- | --- | --- | --- | -------- | --- | --- | --- | --- | --- | --- | --- |
teininterfaceswhilemodulatingminimuminteractions.(a)Thepercentageofcom- The F scores assumed that the true positives were the shared
1
mon interface residues identified by both INTERCAAT and the other database predictedresidues,thefalsepositivesweretheuniqueresiduespre-
dividedbythetotalamountofresiduesidentifiedbytheotherdatabasefordifferent
dictedbyINTERCAATandthefalsenegativesweretheuniqueresi-
minimuminteractioncutoffs.(b)F1scores
|     |     |     |     |     |     | dues predicted | by  | another | database. | As the | minimum | interactions |     |
| --- | --- | --- | --- | --- | --- | -------------- | --- | ------- | --------- | ------ | ------- | ------------ | --- |
requirementwasincreased,weconsistentlyobservedthattherecall
interactionbetweenthequerychainandtheinteractingchain(s),the ofINTERCAATdecreasedandtheprecisionincreased.Aminimum
| distance between | the interacting | atoms | and | the assigned | atom |     |     |     |     |     |     |     |     |
| ---------------- | --------------- | ----- | --- | ------------ | ---- | --- | --- | --- | --- | --- | --- | --- | --- |
interactionrequirementoffouremergedastheoptimalbalancebe-
classes.Thecompatibilitymatrix,ifdisplayed,showseachinterface
|                                                         |     |     |     |     |     | tweenprecisionandrecall.TheresultingF |     |     |     |     | scoresofINTERCAAT |     |     |
| ------------------------------------------------------- | --- | --- | --- | --- | --- | ------------------------------------- | --- | --- | --- | --- | ----------------- | --- | --- |
| residueinthequerychainandthecorrespondingnumberofatomic |     |     |     |     |     |                                       |     |     |     |     | 1                 |     |     |
againsttheotherdatabaseswasapproximately0.8,asexpected(Gil
interactions.
andFiser,2019).
Theatomicclassofanatomisnotpredeterminedbasedonthe
residueitbelongsto.Instead,theprogramdeterminesitsclass-based
| solely on its | coordinates and | particular | atom | type. Therefore, | the | Funding |     |     |     |     |     |     |     |
| ------------- | --------------- | ---------- | ---- | ---------------- | --- | ------- | --- | --- | --- | --- | --- | --- | --- |
programisabletoclassifyatomsfrommostmoleculesincludingpro-
teins,DNA,RNA,etc.INTERCAATcurrentlyrecognizescommon ThisworkwassupportedbyNIHgrantsGM136357andAI141816.
biologicalatoms:C,N,O,P,S,Cl,FandBr.Ifanunknownatomis
inputintoINTERCAATitwillassignitanarbitraryVanderWaal ConflictofInterest:Theauthorsdeclarenoconflictofinterest.
radiusequalto1.8angstromsandassignitsclassas‘?’.Anyatom
| withclass ‘?’ | willbe considereduniversally |     | compatible.It |     | is upto |     |     |     |     |     |     |     |     |
| ------------- | ---------------------------- | --- | ------------- | --- | ------- | --- | --- | --- | --- | --- | --- | --- | --- |
References
| the user to determine | if the interaction |     | makes | sense or update | the |     |     |     |     |     |     |     |     |
| --------------------- | ------------------ | --- | ----- | --------------- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
scripttohandlenewatomtypes. Barber,C.B. et al. (1996) The quickhull algorithm for convex hulls. ACM
INTERCAAT consists of two programs written in python ver- Trans.Math.Softw.,22,469–483.
sion3.8.6,an.iniconfigurationfileandtheqhullpackage.Thetwo Berman,H.M.etal.(2000)TheProteinDataBank.NucleicAcidsRes.,28,
| pythonfilescontainthemainscriptandthefunctions.Theconfigur- |     |     |     |     |     | 235–242. |     |     |     |     |     |     |     |
| ----------------------------------------------------------- | --- | --- | --- | --- | --- | -------- | --- | --- | --- | --- | --- | --- | --- |
ation file must be changed to specify the path to call qhull. Qhull Gil,N.andFiser,A.(2019)Thechoiceofsequencehomologsincludedinmul-
softwarecalculatestheVoronoitessellationbetweenatoms(Barber,
tiplesequencealignmentshasadramaticimpactonevolutionaryconserva-
1996).Iftheuserdoesnothavetheqhullprogramorprefersnotto tionanalysis.Bioinformatics,35,12–19.
useit,itcanbespecifiedintheconfigurationfiletoruntheVoronoi Poupon,A.(2004)VoronoiandVoronoi-relatedtessellationsinstudiesofpro-
calculation using python instead. The downside to this is that the teinstructureandinteraction.Curr.Opin.Struct.Biol.,14,233–241.
programwillrunmuchslower.Forconvenience,bothpythonscripts Richards,F.M.(1974)Theinterpretationofproteinstructures:totalvolume,
arewellcommented. groupvolumedistributionsandpackingdensity.J.Mol.Biol.,82,1–14.
Sobolev,V.etal.(1999)Automatedanalysisofinteratomiccontactsinpro-
teins.Bioinformatics,15,327–332.
3Implementation Vreven,T.etal.(2015)Updatestotheintegratedprotein–proteininteraction
benchmarks:dockingbenchmarkversion5andaffinitybenchmarkversion
Benchmarkingisnotreallypossibleinthesensethatthereisnogold 2.J.Mol.Biol.,427,3031–3041.
standardofinterfacedefinitionsavailable.However,threedifferent
Yang,J.etal.(2013)BioLiP:asemi-manuallycurateddatabaseforbiologically
| databaseswereutilizedtocompareresultsofINTERCAAT. |     |     |     |     | These |          |                |     |               |         |       |       |     |
| ------------------------------------------------- | --- | --- | --- | --- | ----- | -------- | -------------- | --- | ------------- | ------- | ----- | ----- | --- |
|                                                   |     |     |     |     |       | relevant | ligand–protein |     | interactions. | Nucleic | Acids | Res., | 41, |
include320,105and125interfacesdefinedbytheBioLiPdatabase
D1096–D1103.
(Yang et al., 2013), the Nox database (Zhu et al., 2006) and the Zhu,H. et al. (2006) NOXclass: prediction of protein–protein interaction
Dockbdatabase(Vrevenetal.,2015),respectively.Althoughallof types.BMCBioinformatics,7,27.

