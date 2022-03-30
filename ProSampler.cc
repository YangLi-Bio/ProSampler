/****************************************************************************
 *                                                                          *
 *                                 Head files                               *
 *                                                                          *
 ****************************************************************************/

#include <ctype.h> // basic functions in C
#include <math.h> // mathematical functions
#include <cmath> // the library used to calculate the CDF function, may not be used here
#include <algorithm> // basic algorithms regarding stl
#include <deque> // bidirectional queue
#include <time.h> // caount the running time
#include <iomanip> // functions used to , e.g., set precision of numeric values and get memory and time
#include <limits.h> // Definition of infinite values, e.g., -Inf and +Inf
#include <numeric> // include the function to calculate the sum of a vector
#include "Markov.cc" // the codes to generate background sequences using Markov model

using namespace std;


/****************************************************************************
 *                                                                          *
 *                      Macros and complex classes                          *
 *                                                                          *
 ****************************************************************************/

#define TRIVIAL_LEN 6 // length of core parts of motifs
#define DEL_LEN 6 // used to delete redundant motifs
#define OVERLAP 8 // used to weight overlaps between motifs
#define PRECISE 4 // the prcision of PSMs
#define MULTIPLE 100
#define MIN(a, b) ((a)<(b)?(a):(b)) // choose the minimum value from two
#define ALPHA 0.01 // the p-value cutoff (without correction) to select significant k-mers
#define BETA 0.05 // the p-value cutoff (without correction) to select sub-significant k-mers
#define A 15

// the array to save the correspondence between p-values and z-scores
int pval_len = 100; // the length of p-value vector
const double pval_zscore[] = {5.99780701500769, 5.88419335480046, 5.81675774012583, 5.76845830463076, 5.73072886823629, 5.6997262034556, 5.67338878890301, 5.65048045428, 5.63020074085708, 5.61200124417479, 5.59548962865789, 5.58037542028118, 5.56643748350308, 5.5535035351021, 5.54143671278319, 5.53012646857746, 5.51948221566578, 5.50942878574218, 5.49990311062366, 5.49085175210435, 5.48222903230124, 5.47399559729405, 5.46611729879819, 5.45856431288329, 5.45131043784548, 5.44433252920094, 5.43761004084737, 5.4311246493012, 5.4248599435764, 5.41880116739794, 5.4129350034887, 5.40724939194336, 5.401733376418, 5.39637697317127, 5.39117105899695, 5.38610727486682, 5.38117794271074, 5.37637599323976, 5.37169490309786, 5.36712863993063, 5.36267161420255, 5.35831863679091, 5.35406488154408, 5.349905852122, 5.34583735254396, 5.3418554609565, 5.3379565062079, 5.3341370468759, 5.33039385244655, 5.3267238863845, 5.32312429087083, 5.31959237301509, 5.31612559237364, 5.31272154962845, 5.30937797629918, 5.30609272537718, 5.3028637627841, 5.29968915956923, 5.29656708477029, 5.29349579887094, 5.29047364779618, 5.28749905739324, 5.28457052835176, 5.28168663152176, 5.2788460035926, 5.27604734310008, 5.27328940673208, 5.27057100590646, 5.2678910035973, 5.26524831138831, 5.2626418867341, 5.26007073041183, 5.25753388414776, 5.25503042840433, 5.25255948031492, 5.25012019175473, 5.24771174753693, 5.24533336372464, 5.24298428604981, 5.24066378843099, 5.23837117158264, 5.23610576170932, 5.23386690927858, 5.23165398786691, 5.22946639307357, 5.22730354149765, 5.22516486977389, 5.22304983366327, 5.22095790719476, 5.21888858185468, 5.21684136582065, 5.21481578323711, 5.21281137352984, 5.21082769075689, 5.20886430299369, 5.20692079175014, 5.20499675141776, 5.20309178874501, 5.20120552233913, 5.19933758219282, 5.19748760923436, 5.19565525489981, 5.19384018072586, 5.19204205796233, 5.19026056720308, 5.18849539803426, 5.18674624869909, 5.18501282577801, 5.18329484388362, 5.18159202536933, 5.17990410005127, 5.17823080494246, 5.17657188399881, 5.17492708787624, 5.17329617369831, 5.17167890483387, 5.1700750506842, 5.16848438647917, 5.16690669308193, 5.16534175680171, 5.1637893692144, 5.16224932699041, 5.16072143172954, 5.15920548980254, 5.15770131219888, 5.1562087143807, 5.15472751614237, 5.15325754147561, 5.15179861843979, 5.1503505790372, 5.14891325909312, 5.14748649814036, 5.14607013930818, 5.14466402921532, 5.14326801786699, 5.14188195855567, 5.14050570776545, 5.13913912507991, 5.13778207309324, 5.13643441732458, 5.13509602613535, 5.1337667706495, 5.13244652467656, 5.13113516463732, 5.12983256949209, 5.1285386206714, 5.12725320200909, 5.12597619967756, 5.1247075021253, 5.12344700001643, 5.12219458617226, 5.12095015551473, 5.11971360501178, 5.11848483362437, 5.11726374225532, 5.11605023369969, 5.11484421259681, 5.1136455853838, 5.11245426025048, 5.11127014709582, 5.11009315748564, 5.10892320461158, 5.10776020325148, 5.10660406973079, 5.10545472188524, 5.10431207902464, 5.10317606189769, 5.10204659265791, 5.10092359483053, 5.09980699328037, 5.09869671418063, 5.09759268498264, 5.0964948343864, 5.09540309231207, 5.09431738987213, 5.09323765934451, 5.09216383414628, 5.09109584880824, 5.09003363895014, 5.08897714125655, 5.08792629345352, 5.08688103428574, 5.08584130349441, 5.08480704179565, 5.0837781908596, 5.08275469328991, 5.08173649260393, 5.08072353321336, 5.07971576040538, 5.07871312032437, 5.077715559954, 5.07672302709984, 5.07573547037243, 5.07475283917076, 5.07377508366615, 5.07280215478663, 5.07183400420158, 5.0708705843069, 5.06991184821042, 5.06895774971779, 5.06800824331861, 5.067063284173, 5.06612282809844, 5.06518683155695, 5.06425525164256, 5.06332804606914, 5.06240517315846, 5.06148659182857, 5.06057226158243, 5.05966214249685, 5.05875619521167, 5.05785438091917, 5.05695666135378, 5.05606299878197, 5.05517335599245, 5.05428769628648, 5.05340598346854, 5.05252818183713, 5.05165425617576, 5.05078417174423, 5.04991789427002, 5.04905538993992, 5.04819662539186, 5.04734156770684, 5.04649018440117, 5.04564244341873, 5.04479831312357, 5.04395776229247, 5.04312076010787, 5.0422872761508, 5.04145728039402, 5.04063074319532, 5.03980763529093, 5.03898792778911, 5.03817159216384, 5.03735860024863, 5.03654892423054, 5.03574253664423, 5.0349394103662, 5.03413951860909, 5.03334283491618, 5.03254933315592, 5.03175898751662, 5.03097177250125, 5.03018766292229, 5.02940663389678, 5.02862866084138, 5.02785371946758, 5.02708178577698, 5.02631283605668, 5.02554684687478, 5.0247837950759, 5.02402365777687, 5.02326641236244, 5.02251203648114, 5.02176050804114, 5.02101180520626, 5.02026590639204, 5.01952279026181, 5.018782435723, 5.01804482192333, 5.01730992824722, 5.01657773431219, 5.01584821996536, 5.01512136527997, 5.01439715055206, 5.01367555629708, 5.01295656324671, 5.01224015234561, 5.0115263047483, 5.01081500181608, 5.01010622511401, 5.00939995640795, 5.00869617766162, 5.00799487103376, 5.00729601887534, 5.00659960372676, 5.00590560831514, 5.00521401555171, 5.00452480852915, 5.00383797051903, 5.00315348496931, 5.00247133550184, 5.00179150590991, 5.00111398015592, 5.00043874236897, 4.99976577684255, 4.99909506803234, 4.99842660055387, 4.99776035918045, 4.9970963288409, 4.99643449461753, 4.99577484174398, 4.99511735560322, 4.99446202172554, 4.99380882578652, 4.99315775360517, 4.99250879114193, 4.99186192449685, 4.9912171399077, 4.99057442374818, 4.98993376252612, 4.98929514288172, 4.98865855158583, 4.98802397553822, 4.98739140176592, 4.9867608174216, 4.9861322097819, 4.98550556624586, 4.98488087433336, 4.98425812168357, 4.98363729605339, 4.98301838531603, 4.98240137745947, 4.98178626058506, 4.98117302290606, 4.98056165274625, 4.97995213853856, 4.97934446882367, 4.97873863224874, 4.97813461756603, 4.97753241363162, 4.97693200940414, 4.97633339394352, 4.97573655640973, 4.97514148606158, 4.9745481722555, 4.97395660444436, 4.97336677217632, 4.97277866509369, 4.97219227293176, 4.97160758551771, 4.97102459276952, 4.9704432846949, 4.96986365139017, 4.96928568303929, 4.96870936991277, 4.96813470236666, 4.96756167084157, 4.96699026586166, 4.96642047803367, 4.965852298046, 4.96528571666767, 4.9647207247475, 4.96415731321311, 4.96359547307004, 4.96303519540089, 4.96247647136437, 4.96191929219451, 4.96136364919975, 4.96080953376211, 4.96025693733639, 4.9597058514493, 4.95915626769871, 4.95860817775278, 4.95806157334926, 4.95751644629465, 4.95697278846347, 4.95643059179748, 4.95588984830498, 4.95535055006004, 4.95481268920179, 4.95427625793372, 4.95374124852294, 4.95320765329955, 4.95267546465587, 4.95214467504588, 4.95161527698444, 4.9510872630467, 4.95056062586745, 4.95003535814043, 4.94951145261779, 4.94898890210939, 4.9484676994822, 4.94794783765973, 4.94742930962141, 4.94691210840198, 4.94639622709097, 4.94588165883204, 4.94536839682248, 4.94485643431265, 4.94434576460537, 4.94383638105545, 4.94332827706909, 4.94282144610338, 4.94231588166579, 4.94181157731362, 4.9413085266535, 4.9408067233409, 4.94030616107963, 4.93980683362133, 4.93930873476499, 4.93881185835651, 4.93831619828815, 4.93782174849815, 4.9373285029702, 4.93683645573304, 4.93634560085995, 4.93585593246837, 4.93536744471942, 4.93488013181748, 4.93439398800974, 4.93390900758584, 4.93342518487737, 4.93294251425754, 4.9324609901407, 4.93198060698197, 4.93150135927688, 4.93102324156091, 4.93054624840915, 4.93007037443587, 4.92959561429423, 4.9291219626758, 4.92864941431027, 4.92817796396505, 4.92770760644489, 4.9272383365916, 4.92677014928361, 4.92630303943566, 4.92583700199849, 4.92537203195844, 4.92490812433715, 4.92444527419122, 4.92398347661188, 4.92352272672468, 4.92306301968917, 4.92260435069854, 4.92214671497936, 4.92169010779126, 4.92123452442661, 4.92077996021023, 4.92032641049907, 4.91987387068194, 4.91942233617923, 4.91897180244257, 4.91852226495461, 4.91807371922868, 4.91762616080854, 4.91717958526812, 4.91673398821121, 4.91628936527122, 4.9158457121109, 4.91540302442209, 4.91496129792544, 4.91452052837016, 4.91408071153379, 4.91364184322192, 4.91320391926793, 4.91276693553278, 4.91233088790476, 4.91189577229922, 4.91146158465835, 4.91102832095096, 4.91059597717223, 4.91016454934348, 4.90973403351194, 4.90930442575054, 4.90887572215766, 4.90844791885695, 4.90802101199708, 4.90759499775151, 4.90716987231835, 4.90674563192004, 4.90632227280325, 4.90589979123859, 4.90547818352045, 4.90505744596679, 4.90463757491894, 4.9042185667414, 4.90380041782162, 4.90338312456985, 4.90296668341893, 4.9025510908241, 4.90213634326278, 4.90172243723446, 4.90130936926044, 4.90089713588369, 4.90048573366865, 4.90007515920108, 4.89966540908784, 4.89925647995676, 4.89884836845644, 4.89844107125609, 4.89803458504535, 4.89762890653414, 4.89722403245248, 4.89681995955034, 4.89641668459744, 4.89601420438315, 4.89561251571629, 4.89521161542498, 4.89481150035648, 4.89441216737703, 4.89401361337175, 4.89361583524441, 4.89321882991733, 4.89282259433122, 4.89242712544503, 4.89203242023582, 4.89163847569859, 4.89124528884616, 4.890852856709, 4.89046117633516, 4.89007024479003, 4.88968005915629, 4.88929061653373, 4.88890191403912, 4.8885139488061, 4.88812671798503, 4.88774021874284, 4.88735444826296, 4.88696940374512, 4.88658508240529, 4.8862014814755, 4.88581859820376, 4.8854364298539, 4.88505497370547, 4.88467422705363, 4.884294187209, 4.88391485149755, 4.88353621726052, 4.88315828185426, 4.88278104265011, 4.88240449703433, 4.88202864240796, 4.8816534761867, 4.88127899580083, 4.88090519869507, 4.88053208232848, 4.88015964417436, 4.87978788172014, 4.87941679246728, 4.87904637393115, 4.87867662364093, 4.87830753913953, 4.87793911798347, 4.87757135774276, 4.87720425600084, 4.87683781035446, 4.87647201841359, 4.8761068778013, 4.87574238615369, 4.87537854111981, 4.8750153403615, 4.87465278155336, 4.87429086238264, 4.87392958054914, 4.87356893376512, 4.87320891975522, 4.87284953625634, 4.87249078101762, 4.87213265180025, 4.8717751463775, 4.87141826253452, 4.87106199806835, 4.87070635078778, 4.87035131851326, 4.86999689907688, 4.8696430903222, 4.86928989010424, 4.86893729628937, 4.86858530675523, 4.86823391939064, 4.86788313209555, 4.86753294278094, 4.86718334936875, 4.8668343497918, 4.86648594199371, 4.86613812392884, 4.8657908935622, 4.86544424886938, 4.86509818783648, 4.86475270846004, 4.86440780874695, 4.8640634867144, 4.86371974038979, 4.8633765678107, 4.86303396702474, 4.86269193608959, 4.86235047307282, 4.86200957605191, 4.86166924311412, 4.86132947235648, 4.86099026188567, 4.86065160981799, 4.86031351427927, 4.85997597340485, 4.85963898533945, 4.85930254823716, 4.85896666026135, 4.85863131958463, 4.85829652438875, 4.85796227286459, 4.85762856321206, 4.85729539364004, 4.85696276236635, 4.85663066761767, 4.85629910762947, 4.85596808064598, 4.85563758492011, 4.8553076187134, 4.85497818029597, 4.85464926794645, 4.85432087995192, 4.8539930146079, 4.85366567021822, 4.85333884509504, 4.85301253755872, 4.85268674593785, 4.85236146856914, 4.85203670379736, 4.85171244997534, 4.85138870546387, 4.85106546863166, 4.85074273785533, 4.85042051151927, 4.85009878801568, 4.84977756574447, 4.84945684311324, 4.84913661853718, 4.8488168904391, 4.8484976572493, 4.84817891740557, 4.84786066935315, 4.84754291154465, 4.84722564244002, 4.84690886050649, 4.84659256421854, 4.84627675205788, 4.84596142251333, 4.84564657408084, 4.84533220526343, 4.84501831457113, 4.84470490052095, 4.84439196163683, 4.8440794964496, 4.84376750349695, 4.84345598132335, 4.84314492848004, 4.84283434352499, 4.84252422502285, 4.84221457154488, 4.84190538166896, 4.84159665397953, 4.84128838706751, 4.84098057953034, 4.84067322997187, 4.84036633700233, 4.84005989923834, 4.83975391530283, 4.83944838382498, 4.83914330344025, 4.83883867279028, 4.83853449052287, 4.83823075529198, 4.83792746575762, 4.83762462058588, 4.83732221844886, 4.83702025802465, 4.83671873799726, 4.83641765705665, 4.8361170138986, 4.83581680722479, 4.83551703574266, 4.83521769816542, 4.83491879321205, 4.8346203196072, 4.8343222760812, 4.83402466137, 4.83372747421516, 4.83343071336382, 4.83313437756862, 4.83283846558773, 4.83254297618479, 4.83224790812885, 4.83195326019439, 4.83165903116126, 4.83136521981464, 4.83107182494502, 4.83077884534819, 4.83048627982516, 4.83019412718219, 4.8299023862307, 4.82961105578728, 4.82932013467364, 4.82902962171661, 4.82873951574806, 4.82844981560491, 4.82816052012911, 4.82787162816756, 4.82758313857212, 4.8272950501996, 4.82700736191168, 4.82672007257492, 4.82643318106071, 4.82614668624526, 4.82586058700957, 4.8255748822394, 4.82528957082521, 4.82500465166222, 4.82472012365027, 4.8244359856939, 4.82415223670224, 4.82386887558903, 4.8235859012726, 4.8233033126758, 4.82302110872602, 4.82273928835513, 4.8224578504995, 4.82217679409992, 4.82189611810161, 4.82161582145418, 4.82133590311163, 4.8210563620323, 4.82077719717884, 4.82049840751822, 4.82021999202166, 4.81994194966468, 4.81966427942697, 4.81938698029247, 4.81911005124929, 4.81883349128969, 4.81855729941007, 4.81828147461095, 4.81800601589694, 4.81773092227672, 4.81745619276301, 4.81718182637257, 4.81690782212614, 4.81663417904846, 4.81636089616823, 4.81608797251807, 4.81581540713454, 4.81554319905808, 4.815271347333, 4.81499985100748, 4.81472870913352, 4.81445792076695, 4.81418748496737, 4.81391740079816, 4.81364766732645, 4.81337828362311, 4.81310924876271, 4.81284056182352, 4.81257222188747, 4.81230422804016, 4.8120365793708, 4.81176927497224, 4.81150231394091, 4.81123569537681, 4.8109694183835, 4.81070348206809, 4.81043788554119, 4.81017262791692, 4.80990770831289, 4.80964312585017, 4.80937887965326, 4.80911496885011, 4.80885139257206, 4.80858814995385, 4.80832524013361, 4.8080626622528, 4.80780041545624, 4.80753849889205, 4.80727691171168, 4.80701565306985, 4.80675472212457, 4.80649411803706, 4.80623383997183, 4.80597388709659, 4.80571425858223, 4.80545495360286, 4.80519597133575, 4.80493731096131, 4.80467897166311, 4.80442095262784, 4.80416325304527, 4.8039058721083, 4.80364880901286, 4.80339206295799, 4.80313563314573, 4.80287951878117, 4.8026237190724, 4.80236823323053, 4.80211306046963, 4.80185820000674, 4.80160365106187, 4.80134941285793, 4.8010954846208, 4.80084186557924, 4.8005885549649, 4.80033555201233, 4.80008285595891, 4.79983046604492, 4.79957838151342, 4.79932660161034, 4.79907512558439, 4.79882395268708, 4.7985730821727, 4.79832251329831, 4.79807224532372, 4.79782227751147, 4.79757260912684, 4.79732323943782, 4.79707416771508, 4.79682539323199, 4.7965769152646, 4.79632873309161, 4.79608084599435, 4.79583325325681, 4.79558595416559, 4.79533894800988, 4.7950922340815, 4.79484581167481, 4.79459968008677, 4.79435383861689, 4.79410828656721, 4.79386302324231, 4.79361804794931, 4.79337335999781, 4.79312895869991, 4.79288484337021, 4.79264101332575, 4.79239746788607, 4.79215420637313, 4.79191122811132, 4.79166853242748, 4.79142611865084, 4.79118398611304, 4.79094213414812, 4.79070056209248, 4.79045926928488, 4.79021825506647, 4.78997751878072, 4.78973705977342, 4.78949687739272, 4.78925697098904, 4.78901733991513, 4.78877798352602, 4.78853890117901, 4.78830009223368, 4.78806155605185, 4.78782329199762, 4.78758529943728, 4.78734757773938, 4.78711012627466, 4.7868729444161, 4.78663603153884, 4.78639938702021, 4.78616301023973, 4.78592690057907, 4.78569105742205, 4.78545548015464, 4.78522016816494, 4.78498512084319, 4.78475033758172, 4.78451581777497, 4.78428156081948, 4.78404756611389, 4.78381383305887, 4.7835803610572, 4.7833471495137, 4.78311419783522, 4.78288150543068, 4.78264907171099, 4.78241689608912, 4.78218497798, 4.78195331680061, 4.78172191196989, 4.78149076290877, 4.78125986904014, 4.78102922978888, 4.7807988445818, 4.78056871284768, 4.78033883401721, 4.78010920752302, 4.77987983279967, 4.77965070928361, 4.77942183641322, 4.77919321362874, 4.77896484037233, 4.77873671608799, 4.77850884022163, 4.77828121222098, 4.77805383153564, 4.77782669761706, 4.77759980991851, 4.77737316789509, 4.77714677100374, 4.77692061870317, 4.77669471045394, 4.77646904571836, 4.77624362396055, 4.7760184446464, 4.77579350724359, 4.77556881122154, 4.77534435605142, 4.77512014120618, 4.77489616616048, 4.77467243039072, 4.77444893337503, 4.77422567459325, 4.77400265352693, 4.77377986965933, 4.77355732247538, 4.77333501146173, 4.77311293610668, 4.77289109590023, 4.77266949033402, 4.77244811890136, 4.7722269810972, 4.77200607641814, 4.77178540436243, 4.77156496442993, 4.77134475612211, 4.77112477894209, 4.77090503239457, 4.77068551598586, 4.77046622922387, 4.77024717161809, 4.77002834267958, 4.76980974192099, 4.76959136885653, 4.76937322300196, 4.76915530387462, 4.76893761099338, 4.76872014387863, 4.76850290205234, 4.76828588503795, 4.76806909236048, 4.76785252354641, 4.76763617812376, 4.76742005562203, 4.76720415557225, 4.7669884775069, 4.76677302095995, 4.76655778546686, 4.76634277056455, 4.76612797579141, 4.76591340068727, 4.76569904479343, 4.76548490765262, 4.76527098880903, 4.76505728780825, 4.76484380419731, 4.76463053752469, 4.76441748734024, 4.76420465319524, 4.76399203464238, 4.76377963123573, 4.76356744253076, 4.76335546808433, 4.76314370745467, 4.76293216020138, 4.76272082588545, 4.7625097040692, 4.76229879431635, 4.76208809619192, 4.76187760926232, 4.76166733309527, 4.76145726725985, 4.76124741132644, 4.76103776486677, 4.76082832745386, 4.76061909866208, 4.76041007806707, 4.76020126524581, 4.75999265977653, 4.75978426123878, 4.75957606921341, 4.75936808328251, 4.75916030302948, 4.75895272803898, 4.75874535789694, 4.75853819219052, 4.75833123050817, 4.75812447243959, 4.7579179175757, 4.75771156550867, 4.75750541583191, 4.75729946814006, 4.75709372202897, 4.75688817709573, 4.75668283293864, 4.75647768915719, 4.75627274535211, 4.75606800112529, 4.75586345607984, 4.75565910982006, 4.75545496195144, 4.75525101208063, 4.75504725981547, 4.75484370476498, 4.75464034653933, 4.75443718474986, 4.75423421900907, 4.75403144893061, 4.75382887412927, 4.75362649422101, 4.7534243088229};
  
  typedef struct str_int_int
  {
    string header;
    char strand;
    int begin;
    int end;
  } st; // a motif site
  
  typedef vector<int> int_vector;
  typedef vector<float> float_vector;
  typedef vector<char> char_vector;
  
  typedef struct kmer_occurr
  {
    vector<string> kmer;
    int_vector occurr;
    float_vector cover1;
    float_vector cover2;
    float_vector score;
    vector<int_vector> site;
  } kmer_set; // contain the k-mers, occurrences, and sites
  
  typedef struct kmer_lmer
  {
    string kmer;
    int occurr;
    float cover1;
    float cover2;
    float score;
    vector<int> site;
    vector<float_vector> left;
    vector<float_vector> right;
  } kmer_all; // k-mers and their flanking strings
  
  typedef struct psm_occurr
  {
    vector<int> set;
    vector<float_vector> mat;
    int occurr;
    string cons;
    float score;
    float cover1;
    float cover2;
  } psm; // PSMs
  
  typedef struct left_psm_right
  {
    vector<int> set;
    vector<float_vector> mat;
    int occurr;
    vector<float_vector> left;
    vector<float_vector> right;
    vector<int> site;
    int begin;
    int end;
    float cover1;
    float cover2;
    float score;
    double pvalue; // chi square p-value
    double padj; // adjusted p-value
  } pre_mtf; // preliminary motifs
  
  typedef struct set_mat_alength_nsites_cons
  {
    vector<int> set;
    vector<float_vector> mat;
    int alength;
    int nsites;
    string cons;
    string rev_cons;
    string deg;
    string rev_deg;
    float score;
    vector<st> site;
    float cover1;
    float cover2;
    double pvalue; // chi square p-value
    double padj; // chi square adjusted p-value
  } mtf; // motifs
  
  typedef struct mat_vec_mat
  {
    vector<float_vector> pfm;
    vector<float> ic;
    vector<float_vector> pssm;
  } spic; // the motif in SPIC formats used to calculate pairwise motif similarity
  
  int f_in_id = -1; // foreground file
  // int f_bg_id = -1; // background file
  int f_out_id = -1; // output file
  int top_num; // the number of preliminary motifs to consider
  int correction_cutoff = 100000; // Cutoff 
  
  int kmer_length = 9; // default k-mer length
  int lmer_length = 9; // the length of flanking regions of each k-mer
  int num_deg = 1; // the number of degenerate positions in a k-mer
  int hd_thr = 1; // the Hamming distance cutoff to merge k-mers into PSMs
  int redundant_thr = 2; // the cutoff to weight the similarity between motifs
  int num_mtf = -1; // the number of motifs to output
  int str_flag = 2; // the number of strands to consider: 1 for positive strands, and 2 for both
  int help_flag = 0; // whether print the helping message
  int num_iter = 200; // the number of iterations in Gibbs sampling
  float sw_thr = 1.8; // the cutoff of edge weights to build the PSM similarity graph
  
  float sum_exp = 0;
  float sum_bg = 0;
  float padj_cutoff = 0.05; // cutoff of adjusted p-values (Bonferroni correction)
  int correction_method = 1; // Metho to perform p-value correction, Bonferroni by default
  
  float thr1 = 0.01; // the adjusted p-value cutoff to select significant k-mers
  float thr2 = 1.96; // the z-score cutoff to determine the motif length
  float thr3 = 0.05; // the adjusted p-value cutoff to select sub-significant k-mers
  
  char alphabet[4] = {'A', 'C', 'G', 'T'};
  map<char, int> r_alphabet;
  float nt_freq[4]; // nucleotide frequencies
  char wild_card[14], reverse_wild[14], begin_pkg[64], end_pkg[64]; // wild cards and recorded time points
  time_t begin_t, end_t; // the beginning and ending time points
  
  
  /****************************************************************************
   *                                                                          *
   *                                 Functions                                *
   *                                                                          *
   ****************************************************************************/
  
  // Free the memory of array
  template <typename T>
  void free_pt(T* tmp_pt)
  {
    if (tmp_pt != NULL)
    {
      delete [] tmp_pt;
      tmp_pt = NULL;
    }
  }
  
  
  // Reverse the order of ranked values in an int array
  void reverse_sort(int id[], int length)
  {
    int *temp_id;
    temp_id = new int [length];
    int i;
    
    for (i = 0; i < length; i ++)
    {
      temp_id[i] = id[length - 1 - i];
    }
    for (i = 0; i < length; i ++)
    {
      id[i] = temp_id[i];
    }
    
    free_pt(temp_id);
  }
  
  
  
  // Print the help message
  void usage()
  {
    cout << "\n=================================================================================================================================================================\n";
    cout << "Usage: ./ProSampler [options]\n";
    cout << "\nParameters:\n";
    cout << "-i: Name of the input file in FASTA format\n";
    // cout << "-b: Name of the background file in FASTA format or order of the Markov " 
    //      << "model to generate background sequences (default: 3; 3rd order Markov model)\n";
    cout << "-d: The cutoff of Hamming Distance between any two k-mers in a PWM (default: 1)\n";
    cout << "-o: Prefix of the names of output files\n";
    cout << "-m: Number of motifs to be output (default: All)\n";
    cout << "-f: Number of cycles of Gibbs Sampling to find each preliminary motif (default: 200)\n";
    cout << "-k: Length of preliminary motifs (default: 9)\n";
    cout << "-l: Length of the flanking l-mers (default: 9)\n";
    cout << "-c: Cutoff of Hamming distance to merge similar k-mers (default: 1)\n";
    cout << "-r: Cutoff of Hamming distance to delete redundant motifs basedn on consensus (default: 1)\n";
    cout << "-p: Number (1 or 2) of strands to be considered (default: 2)\n";
    cout << "-P: cutoff of adjusted p-value for motifs (default: 0.05)\n";
    cout << "-t: Cutoff of p-value to choose significant k-mers (default: 0.01)\n";
    cout << "-w: Cutoff of p-value to choose sub-significant k-mers (default:0.05)\n";
    cout << "-z: Cutoff of p-value to extend preliniary motifs(default: 1.96)\n";
    cout << "-s: Cutoff of SW score to construct graph (default: 1.80)\n";
    cout << "-C: Method to perform p-value correction (default: 1 - Bonferroni; 2 - Benjamini-Hochberg)\n";
    cout << "-R: Cutoff of data size to perform p-value correction (100000 by default)\n";
    cout << "-h: Print this message (default: 0)\n";
    cout << "=================================================================================================================================================================\n";
  }
  
  
  // Parse the input options
  void parse_opt(int argc, char** argv)
  {
    int i;
    
    if (argc < 3)
    {
      cout << "Error: The number of input parameters is wrong. Please check the input parameters again.\n";
      usage();
      exit(0);
    }
    
    for (i = 1; i < argc - 1; i ++)
    {
      if (strcmp(argv[i], "-i") == 0) {f_in_id = i + 1; i ++;}
      // else if (strcmp(argv[i], "-b") == 0) {f_bg_id = i + 1; i ++;}
      else if (strcmp(argv[i], "-o") == 0) {f_out_id = i + 1; i ++;}
      else if (strcmp(argv[i], "-d") == 0) {num_deg = atoi(argv[i + 1]); i ++;}
      else if (strcmp(argv[i], "-k") == 0) {kmer_length = atoi(argv[i + 1]); i ++;}
      else if (strcmp(argv[i], "-l") == 0) {lmer_length = atoi(argv[i + 1]); i ++;}
      else if (strcmp(argv[i], "-m") == 0) {num_mtf = atoi(argv[i + 1]); i ++;}
      else if (strcmp(argv[i], "-c") == 0) {hd_thr = atoi(argv[i + 1]); i ++;}
      else if (strcmp(argv[i], "-r") == 0) {redundant_thr = atoi(argv[i + 1]); i ++;}
      else if (strcmp(argv[i], "-p") == 0) {str_flag = atoi(argv[i + 1]); i ++;}
      else if (strcmp(argv[i], "-f") == 0) {num_iter = atoi(argv[i + 1]); i ++;}
      else if (strcmp(argv[i], "-t") == 0) {thr1 = atof(argv[i + 1]); i ++;}
      else if (strcmp(argv[i], "-w") == 0) {thr3 = atof(argv[i + 1]); i ++;}
      else if (strcmp(argv[i], "-z") == 0) {thr2 = atof(argv[i + 1]); i ++;}
      else if (strcmp(argv[i], "-h") == 0) {help_flag = atoi(argv[i + 1]); i ++;}
      else if (strcmp(argv[i], "-R") == 0) {correction_cutoff = atoi(argv[i + 1]); i ++;}
      else if (strcmp(argv[i], "-s") == 0) {sw_thr = atof(argv[i + 1]); i ++;}
      else if (strcmp(argv[i], "-P") == 0) {padj_cutoff = atof(argv[i + 1]); i ++;}
      else if (strcmp(argv[i], "-C") == 0) {correction_method = atoi(argv[i + 1]); i ++;}
    }
    
    cout << "Finished parsing input parameters to the program.\n";
  }
  
  
  // Unify all the nucleotides to upper cases
  void de_lower(sequence& seq)
  {
    int i, j;
    
    for (i = 0; i < seq.strand.size(); i ++)
    {
      for (j = 0; j < seq.strand[i].length(); j ++)
      {
        if (islower(seq.strand[i][j]))
        {
          seq.strand[i][j] = toupper(seq.strand[i][j]); // transform lower case to UPPER case
        }
      }
    }
  }
  
  
  // Substitute all the unknown letters to "N"
  void de_other(sequence& seq)
  {
    for (int i = 0; i < seq.strand.size(); i ++)
    {
      for (int j = 0; j < seq.strand[i].length(); j ++)
      {
        if (seq.strand[i][j] != 'A' && seq.strand[i][j] != 'C' 
              && seq.strand[i][j] != 'G' && 
                seq.strand[i][j] != 'T')
        {
          seq.strand[i][j] = 'N'; // all nonregular characters are transformed to 'N'
        }
      }
    }
    
    cout << "Finished getting rid of the unknown letters in FASTA file.\n";
  }
  
  
  // Concatenate two strings
  void copy_str(string& small, string& big)
  {
    for (int i = 0; i < small.length(); i ++)
    {
      if (small[i] != ' ')
      {
        big += small[i]; 
      }
    }
  }
  
  
  // Get the reversed strand for a sequence
  void seq_reverse(sequence& seq, vector<string>& strand)
  {
    int i, j;
    string str;
    strand.clear();
    
    for (i = 0; i < seq.strand.size(); i ++)
    {
      str.assign(seq.strand[i].rbegin(), seq.strand[i].rend());
      for (j = 0; j < str.length(); j ++)
      {
        switch (str[j]) // substitute each nucleotide with its supplement
        {
        case 'A':
          str[j] = 'T';
          break;
        case 'C':
          str[j] = 'G';
          break;
        case 'G':
          str[j] = 'C';
          break;
        case 'T':
          str[j] = 'A';
          break;
        default:
          break;
        }
      }
      strand.push_back(str);
    }
  }
  
  
  // Generate the final version of sequences for motif finding
  void seq_generate(sequence& seq_pos, sequence& seq_final, 
                    int str_flag)
  {
    if (str_flag == 1)
    {
      seq_final = seq_pos; // keep the same
    }
    else if (str_flag == 2)
    {
      vector<string> strand;
      seq_reverse(seq_pos, strand);
      
      if (seq_pos.strand.size() != strand.size()) // the number of sequences in positive and negative
        // directions are not equal
      {
        
        cerr << "Error: The number of positive strands of " 
             << "the sequences is not equal to the number of negative strands.\n";
        exit(1);
      }
      seq_final.name = seq_pos.name;
      seq_final.strand.clear();
      
      for (int i = 0; i < seq_final.name.size(); i ++)
      {
        seq_final.strand.push_back(seq_pos.strand[i]);
        seq_final.strand.push_back(strand[i]);
      }
      
      if (seq_final.name.size() * 2 != seq_final.strand.size())
      {
        cerr << "Error: The number of sequence names " 
             << "mismatches that of sequence strands.\n";
        exit(1);
      }
    }
    else
    {
      cerr << "Error: The parameter str_flag should be chosen from 1 or 2.\n" 
           << "Please check the parameters again.\n";
      exit(1);
    }
  }
  
  
  // Load data from a file into a sequence object
  void load_data(char* file_name, sequence& seq)
  {
    seq.name.clear();
    seq.strand.clear();
    
    ifstream file_op(file_name); // open a file handle to read sequences
    if (!file_op) // failed in opening the file
    {
      cerr << "Error: Can't open file " << file_name << "!\n";
      exit(1);
    }
    
    string line; // one line in the file handle
    string merge_ln; // one sequence saved in multiple lines will be merged into one
    
    while (!file_op.eof()) // while the files has not ended
    {
      getline(file_op, line); // read one line
      
      if (line[0] == '>') // '>' represents a line corresponding to the sequence name
      {
        if (!merge_ln.empty()) // if the strand line is not empty
        {
          seq.strand.push_back(merge_ln); // save the merged line for the previous sequence
        }
        
        seq.name.push_back(line.substr(1)); // save the sequence name
        merge_ln.clear(); // launch a new sequence
      }
      else if (!line.empty()) // if this is a nonempty strand line
      {
        copy_str(line, merge_ln); // concatenate 
        // concatenate one string at the end of another string
      }
      else // if this is an empty line, we will skip it.
      {
        continue;
      }
    }
    
    if (!merge_ln.empty())
    {
      seq.strand.push_back(merge_ln);
    }
    
    if (seq.name.size() != seq.strand.size()) // if the number of sequence names is unequal with 
      // that of sequence strands
    {
      cerr << "Error: The number of sequence name lines is not " 
           << "equal to the number of sequence lines!\n";
      exit(1);
    }
    
    cout << "Finished loading FASTA file " << file_name 
         << " into sequence variable.\n";
    cout << "There're altogether " << seq.name.size() 
         << " FASTA sequences inputted into this program.\n";
  }
  
  
  // Free the space of sequences
  void free_seq(sequence& seq)
  {
    for (int i = 0; i < seq.name.size(); i ++)
    {
      string().swap(seq.name[i]); // free the space to save the sequences
    }
    
    vector<string>().swap(seq.name);
    
    for (int i = 0; i < seq.strand.size(); i ++)
    {
      string().swap(seq.strand[i]);
    }
    
    vector<string>().swap(seq.strand);
  }
  
  
  // Calculate the frequencies of nucleotides
  void nt_stat(sequence& seq_nondeg, float* nt_freq)
  {
    int i, j;
    float sum = 0;
    
    for (i = 0; i < 4; i ++)
    {
      nt_freq[i] = 0;
    }
    
    for (i = 0; i < seq_nondeg.strand.size(); i ++)
    {
      for (j = 0; j < seq_nondeg.strand[i].length(); j ++)
      {
        switch (seq_nondeg.strand[i][j])
        {
        case 'A':
          nt_freq[0]++;
          break;
        case 'C':
          nt_freq[1]++;
          break;
        case 'G':
          nt_freq[2]++;
          break;
        case 'T':
          nt_freq[3]++;
          break;
        default:
          break;
        }
      }
    }
    
    for (i = 0; i < 4; i ++)
    {
      sum += nt_freq[i];
    }
    
    for (i = 0; i < 4; i ++)
    {
      nt_freq[i] /= sum;
    }
    
    cout << "Finished counting the nucleotide frequencies of the input FASTA file.\n";
  }
  
  
  // Count the occurrence numbers of k-mers
  void kmer_count(map<string, int_vector>& kmer_site, 
                  sequence& seq_final, int& lmer_length)
  {
    int i, j;
    
    kmer_site.clear();
    
    if (str_flag == 1)
    {
      for (i = 0; i < seq_final.strand.size(); i ++)
      {
        if (seq_final.strand[i].length() < 
          kmer_length + 2 * lmer_length)
        {
          continue;
        }
        
        for (j = lmer_length; j <= seq_final.strand[i].
        length() - kmer_length - lmer_length; j ++)
        {
          if (seq_final.strand[i].substr(j, kmer_length).
                find('N') != seq_final
                .strand[i].substr(j, kmer_length).npos)
          {
            continue;
          }
          
          kmer_site[seq_final.strand[i].substr(j, kmer_length)]
                             .push_back(i);
          kmer_site[seq_final.strand[i].substr(j, kmer_length)]
                             .push_back(j);
        }
      }
    }
    else
    {
      for (i = 0; i < seq_final.strand.size(); i += 2)
      {
        if (seq_final.strand[i].length() < kmer_length
              + 2 * lmer_length)
        {
          continue;
        }
        
        for (j = lmer_length; j <= seq_final.strand[i].
        length() - kmer_length - lmer_length; j ++)
        {
          if (seq_final.strand[i].substr(j, kmer_length)
                       .find('N') != seq_final.strand[i]
                       .substr(j, kmer_length).npos)
          {
            continue;
          }
          
          kmer_site[seq_final.strand[i]
                             .substr(j, kmer_length)].push_back(i);
          kmer_site[seq_final.strand[i]
                             .substr(j, kmer_length)].push_back(j);
        }
      }
    }
    
  }
  
  
  // Get rid of the k-mers with low complexity (length: 5)
  int if_same(string str)
  {
    int i, flag;
    flag=0;
    
    for (i = 1; i < str.length(); i ++)
    {
      if (str[0] != str[i])
      {
        flag = 1;
        break;
      }
    }
    
    return flag;
  }
  
  
  // Check complexity of a string
  int check_complex(string str)
  {
    int flag = 0;
    
    map<char, int> map_str;
    
    for (int i = 0; i < str.length(); i ++)
    {
      map_str[str[i]] ++;
    }
    
    if (map_str.size() >= 2)
    {
      flag = 1;
    }
    
    return flag;
  }
  
  
  // Get rid of the k-mers with low complexity (unfixed length)
  int de_trivial(string str)
  {
    int flag = 1;
    
    for (int i = 0; i <= str.length() - TRIVIAL_LEN; i ++)
    {
      if (!if_same(str.substr(i, TRIVIAL_LEN)))
      {
        flag = 0;
        break;
      }
    }
    
    if(flag == 1)
    {
      flag = check_complex(str);
      
    }
    
    return flag;
  }
  
  
  // Get rid of the k-mers with low complexity (unfixed length)
  int last_trivial(string str)
  {
    int flag = 1;
    
    for (int i = 0; i <= str.length() - DEL_LEN; i ++)
    {
      if (!if_same(str.substr(i, DEL_LEN)))
      {
        flag = 0;
        break;
      }
    }
    
    if (flag == 1)
    {
      flag = check_complex(str);
      
    }
    
    return flag;
  }
  
  
  // Perform two proportion z-test
  float ztest(float& cover1, float& cover2, float& sum1, 
              float& sum2)
  {
    float p1 = cover1/sum1;
    float p2 = cover2/sum2;
    float p = (cover1 + cover2) / (sum1 + sum2);
    float up = p1 - p2;
    float down = sqrt(p * (1 - p) * (1 / sum1 + 1 / sum2));
    float res = up / down;
    return res;
  }
  
  
  // Initialize the array
  void init_array(int* flag_ar, int len, int start)
  {
    for(int i = 0; i < len; i ++)
    {
      flag_ar[i] = start;
    }
  }
  
  
  // Transform site information into coverage information
  float site2cover(vector<int>& site, int seq_num)
  {
    int* seq_ar;
    seq_ar = new int [seq_num];
    init_array(seq_ar, seq_num, 0);
    float sum = 0;
    
    if (str_flag == 2)
    {
      for (int i = 0; i < site.size(); i += 2)
      {
        seq_ar[site[i] / 2] = 1;
      }
    }
    else
    {
      for (int i = 0; i < site.size(); i += 2)
      {
        seq_ar[site[i]] = 1;
      }
    }
    
    for (int i = 0; i < seq_num; i ++)
    {
      sum += seq_ar[i];
    }
    
    free_pt(seq_ar);
    
    return sum;
  }
  
  
  // Calculate the CDF 
  double normalCDF(double x)
  {
    //int n = (int)(x / 1e-08) // convert p-values into intervals of integers
    
    int n = (int)(x / 1e-08);
    double zscore;
    
    if (n < 0)
    {
      zscore = pval_zscore[0];
    }
    else if (n >= 100)
    {
      zscore = pval_zscore[99];
    }
    else
    {
      zscore = pval_zscore[(int)(x / 1e-08)];
    }
    
    return(zscore);
  }
  
  
  // Convert z-scores to p-values
  float zscore_to_pval(float zscore) {
    
    // zscore : the z-score
    
    int id = 0;
    for (int i = 0; i < pval_len; i ++) {
      if (zscore <= pval_zscore[i]) {
        id = i;
        break;
      }
    }
    
    return(id * 1e-08);
  }
  
  
  // BH correction
  void BH_correct(kmer_set& major_set, kmer_set& minor_set, float thr1, float thr3) {
    
    // major_set : all the k-mers
    // minor_set : empty set
    // thr1 : the p-value cutoff to select significant k-mers
    // thr3 : the p-value cutoff to select subsignificant k-mers
    
    // Select significant k-mers
    float cutoff = thr1 / major_set.kmer.size();
    int id = 0;
    for (int i = 0; i < major_set.kmer.size(); i ++) {
      float pval = zscore_to_pval(major_set.score[i]); // convert z-score to the p-value
      if (pval > cutoff * (i + 1)) {
        id = i - 1;
        break;
      }
    }
    if (id < 0) {
      id = 0;
    }
    kmer_set tmp_set; // a temporary set of k-mers
    // tmp_set.reserve(id + 1);
    for (int i = 0; i <= id; i ++) {
      tmp_set.kmer.push_back(major_set.kmer[i]);
      tmp_set.occurr.push_back(major_set.occurr[i]);
      tmp_set.cover1.push_back(major_set.cover1[i]);
      tmp_set.cover2.push_back(major_set.cover2[i]);
      tmp_set.score.push_back(major_set.score[i]);
      tmp_set.site.push_back(major_set.site[i]);
    }
    major_set = tmp_set;
    
    
    // Select subsignificant k-mers
    cutoff = thr3 / major_set.kmer.size();
    int id2 = id + 1;
    for (int i = id2; i < major_set.kmer.size(); i ++) {
      float pval = zscore_to_pval(major_set.score[i]); // convert z-score to the p-value
      // cout << pval << "\n";
      if (pval > cutoff * (i + 1)) {
        id2 = i - 1;
        break;
      }
    }
    if (id2 < 0) {
      id2 = 0;
    }
    kmer_set tmp_set2; // a temporary set of k-mers
    // tmp_set2.reserve(id2 - id + 1);
    for (int i = id + 1; i <= id2; i ++) {
      tmp_set2.kmer.push_back(major_set.kmer[i]);
      tmp_set2.occurr.push_back(major_set.occurr[i]);
      tmp_set2.cover1.push_back(major_set.cover1[i]);
      tmp_set2.cover2.push_back(major_set.cover2[i]);
      tmp_set2.score.push_back(major_set.score[i]);
      tmp_set2.site.push_back(major_set.site[i]);
    }
    minor_set = tmp_set2;
  }
  
  
  // Select significant k-mers
  void choose_kmer(map<string, int_vector>& map_exp, 
                   map<string, float>& map_bg, kmer_set& major_set, int correction_method, 
                   int seq_num, map<string, int>& fa_flag, kmer_set& minor_set, 
                   float thr1, float thr3)
  {
    // map_exp : the hash table of k-mers in the positive dataset
    // correction_method : method for correction: 0 - no correction; 1 - Bonforoni; 2 - BH
    // map_bg : the k-mer frequencies in background
    // major_set : significant k-mers
    // seq_num : the number of sequences
    // fa_flag : flag for k-mers
    // minor_set : the set of subsignificant k-mers
    // thr1 : adjusted p-value cutoff for significant k-mers, 0.01 by default
    // thr3 : adjusted p-value cutoff  for subsignificant k-mers, 0.05 by default
    
    major_set.kmer.clear();
    major_set.cover1.clear();
    major_set.cover2.clear();
    major_set.site.clear();
    major_set.score.clear();
    
    float tmp_cover1, tmp_cover2;
    float z_sc;
    sum_exp = 0;
    sum_bg = 0;
    
    if (correction_method == 1) { // Bonferoni correction 
      thr1 /= double(map_exp.size());
      thr3 /= double(map_exp.size());
      thr1 = normalCDF(thr1); // find the z-score cutoff to select significant k-mers
      thr3 = normalCDF(thr3); // find the z-score cutoff to select significant k-mers
    }
    
    cout << "The z-score cutoff to choose significant k-mers is " << thr1 << "\n";
    cout << "The z-score cutoff to choose significant k-mers is " << thr3 << "\n";
    
    for (map<string, int_vector>::iterator i = 
         map_exp.begin(); i != map_exp.end(); i ++)
    {
      if (fa_flag[i -> first] == 0)
      {
        continue;
      }
      
      sum_exp += (i -> second).size() / 2;
    }
    
    for (map<string, float>::iterator i = map_bg.begin(); 
         i != map_bg.end(); i ++)
    {
      sum_bg += int(i -> second);
    }
    
    for (map<string, int_vector>::iterator i = map_exp.begin(); 
         i != map_exp.end(); i ++)
    {
      if (fa_flag[i -> first] == 0)
      {
        continue;
      }
      
      tmp_cover1 = (i -> second).size() / 2;
      tmp_cover2 = int(map_bg[i -> first]);
      
      // cout << i -> first << " " << tmp_cover1 << " " << tmp_cover2 << "\n";
      // exit(0);
      
      z_sc = ztest(tmp_cover1, tmp_cover2, sum_exp, sum_bg); // two proportion z-test
      cout << z_sc << "\n";

      // Not BH correction
      if (correction_method != 2) {
        if (z_sc > thr1)
        {
          major_set.kmer.push_back(i -> first);
          major_set.occurr.push_back(tmp_cover1);
          major_set.cover1.push_back(tmp_cover1);
          major_set.cover2.push_back(tmp_cover2);
          major_set.site.push_back(i -> second);
          major_set.score.push_back(z_sc);
          
          // cout << i -> first << " " << tmp_cover1 << " " << tmp_cover2 << "\n";
          // exit(0);
        }
        else if (z_sc > thr3)
        {
          minor_set.kmer.push_back(i -> first);
          minor_set.occurr.push_back(tmp_cover1);
          minor_set.cover1.push_back(tmp_cover1);
          minor_set.cover2.push_back(tmp_cover2);
          minor_set.site.push_back(i -> second);
          minor_set.score.push_back(z_sc);
          
          // cout << i -> first << " " << tmp_cover1 << " " << tmp_cover2 << "\n";
          // exit(0);
        }
      } else {
        major_set.kmer.push_back(i -> first);
        major_set.occurr.push_back(tmp_cover1);
        major_set.cover1.push_back(tmp_cover1);
        major_set.cover2.push_back(tmp_cover2);
        major_set.site.push_back(i -> second);
        major_set.score.push_back(z_sc);
      }
    }
    
    // cout << major_set.kmer[0] << " " << major_set.cover1[0] << " " << major_set.cover2[0] << "\n";
    // exit(0);
    
    cout << "there are altogether " << major_set.kmer.size() << 
      " significant " << kmer_length << "-mers left.\n";
    cout << "there are altogether " << minor_set.kmer.size() << 
      " sub-significant " << kmer_length << "-mers left.\n";
  }
  
  
  // Get the reverse strand of a string
  string get_rev(string str)
  {
    string rev_str;	// The reverse string
    
    for (int i = str.length() - 1; i >= 0; i --)
    {
      rev_str += alphabet[3 - r_alphabet[str[i]]];
    }
    
    return rev_str;
  }
  
  
  // Check whether the string is palindrome
  int check_parlindrome(string str)
  {
    string rev_str = get_rev(str);	// The reverse string
    
    int flag = 0;
    
    if (str == rev_str)
    {
      flag = 1;
    }
    
    return flag;
  }
  
  
  // Reverse the sites of reverse supplementary strands
  void reverse_st(vector<int>& positive, vector<int>& negative, 
                  vector<string>& strand)
  {
    negative = positive;
    
    for (int i = 0; i < negative.size() - 1; i += 2)
    {
      negative[i]++;
    }
    
    for (int i=1; i<negative.size(); i += 2)
    {
      negative[i]=strand[negative[i-1]].length()-negative[i]-kmer_length;
    }
  }
  
  
  // Combine reverse complementary k-mers
  void comb_kmer(map<string, int_vector>& kmer_site, map<string, int>& fa_flag, vector<string>& positive, 
                 map<string, float>& kmer_bg, map<string, int>& bg_flag)
  {
    // kmer_site : the hash table of k-mers in positive data
    // fa_flag : to indicate whether we have counted this k-mer in the positive dataset
    // positive : the positive data composed of all sequences
    // kmer_bg : the frequencies of k-mers
    // bg_flag : to indicate whether we have counted this k-mer in the negative dataset
    
    for(map<string, int_vector>::iterator it = kmer_site.begin(); it != kmer_site.end(); it ++)
    {
      fa_flag[it -> first] = 1; // We have counted this k-mer
    }
    
    for(map<string, int_vector>::iterator it = kmer_site.begin(); it != kmer_site.end(); it ++)
    {
      if(fa_flag[it -> first] == 0 || check_parlindrome(it -> first))
      {
        continue;
      }
      
      string rev_str = get_rev(it -> first);
      
      if(kmer_site.count(rev_str) > 0)
      {
        fa_flag[rev_str] = 0;
        vector<int> rev_site;
        rev_site.reserve(kmer_site[rev_str].size());
        reverse_st(kmer_site[rev_str], rev_site, positive);
        (it->second).insert((it->second).end(), rev_site.begin(), rev_site.end());
        kmer_site[rev_str].clear();
      }
    }
    
    bg_flag = fa_flag;
    
    for(map<string, int>::iterator it = bg_flag.begin(); it != bg_flag.end(); it ++)
    {
      if(it -> second == 0)
      {
        continue;
      }
      
      string rev_str = get_rev(it->first);
      
      if(kmer_bg.count(rev_str) > 0) // Add the reserved string
      {
        // vector<int> rev_site;
        // rev_site.reserve(kmer_bg[rev_str].size());
        // reverse_st(kmer_bg[rev_str], rev_site, negative);
        // kmer_bg[it -> first].insert(kmer_bg[it -> first].end(), rev_site.begin(), rev_site.end());
        kmer_bg[it -> first] += kmer_bg[rev_str]; 
        // kmer_bg[rev_str].clear();
      }
    }
    
    cout << "There are " << kmer_site.size() << " " << kmer_length << 
      "-mers after combining reverse-complementary ones.\n";
  }
  
  
  // Sort the chosen k-mers
  void kmer_sort(kmer_set& chosen_kmers)
  {
    float *kmer_score;
    int *kmer_id1;
    int num=chosen_kmers.score.size();
    kmer_score=new float [num];
    kmer_id1=new int [num];
    int i;
    
    for(i=0; i<num; i++)
    {
      kmer_score[i]=chosen_kmers.score[i];
      kmer_id1[i]=i;
    }
    
    quickSort(kmer_score, 0, num-1, kmer_id1);
    reverse_sort(kmer_id1, num);
    
    kmer_set temp_kmers;
    for(i=0; i<num; i++)
    {
      temp_kmers.kmer.push_back(chosen_kmers.kmer[kmer_id1[i]]);
      temp_kmers.score.push_back(chosen_kmers.score[kmer_id1[i]]);
      temp_kmers.cover1.push_back(chosen_kmers.cover1[kmer_id1[i]]);
      temp_kmers.cover2.push_back(chosen_kmers.cover2[kmer_id1[i]]);
      temp_kmers.occurr.push_back(chosen_kmers.occurr[kmer_id1[i]]);
      temp_kmers.site.push_back(chosen_kmers.site[kmer_id1[i]]);
    }
    
    chosen_kmers=temp_kmers;
    
    delete [] kmer_id1;
    delete [] kmer_score;
    
    cout<<"Finished sorting the significant "<<kmer_length<<"-mers based on z-scores.\n";
  }
  
  
  // Sort all the PSMs in decreasing order of PSM scores and simultaneously delete the trivial PSMs
  void psm_sort(vector<psm>& psm_set)
  {
    int *psm_id1;
    float *psm_scores;
    int num=psm_set.size();
    psm_scores=new float [num];
    psm_id1=new int [num];
    
    int i;
    
    for(i=0; i<num; i++)
    {
      psm_scores[i]=psm_set[i].score;
      psm_id1[i]=i;
    }
    
    quickSort(psm_scores, 0, num-1, psm_id1);
    reverse_sort(psm_id1, num);
    
    {
      vector<psm> tmp_psm;
      for(i=0; i<num; i++)
      {
        tmp_psm.push_back(psm_set[psm_id1[i]]);
      }
      
      psm_set.swap(tmp_psm);
    }
    
    delete [] psm_id1;
    delete [] psm_scores;
    
    cout<<"Finished sorting the PSMs in PSM set.\n";
    cout<<"Finished kicking out the trivial PSMs.\n"<<"There're "<<psm_set.size()<<" PSMs left."<<endl;
  }
  
  
  // Initialize a flag vector
  void generate_flag(int length, vector<int>& flag)
  {
    int i;
    
    flag.reserve(length);
    
    for(i=0; i<length; i++)
    {
      flag.push_back(0);
    }
  }
  
  
  // Mask_kmer used to mask the selected k-mers
  void mask_kmer(vector<int>& set, vector<int>& kmer_flag)
  {
    int i;
    
    for(i=0; i<set.size(); i++)
    {
      kmer_flag[set[i]]=1;
    }
  }
  
  
  // Find the id of next seed for partitioning k-mers
  int find_seed(vector<int>& kmer_flag)
  {
    int out_id=-1;
    int i;
    
    for(i=0; i<kmer_flag.size(); i++)
    {
      if(kmer_flag[i]==0)
      {
        out_id=i;
        break;
      }
    }
    
    return out_id;
  }
  
  
  // Get the Hamming Distance between two strings
  int comp_hd(string str1, string str2)
  {
    if(str1.length()!=str2.length())
    {
      cerr<<"Error: only two strings with the same length can be compared to get Hamming Distance.\n";
      exit(0);
    }
    
    int hd=0;
    
    for(int i=0; i<str1.length(); i++)
    {
      if(str1[i]!=str2[i])
      {
        hd++;
      }
    }
    
    return hd;
  }
  
  
  // Transform all significant k-mers into PSMs
  void kmer2psm(vector<kmer_all>& kmer_final, vector<psm>& psm_set, int& top_num)
  {
    psm_set.clear();
    psm tmp_psm;
    
    for(int i=0; i<top_num; i++)
    {
      tmp_psm.set.clear();
      tmp_psm.set.push_back(i);
      psm_set.push_back(tmp_psm);
    }
    
    for(int i=1; i<top_num; i++)
    {
      for(int j=0; j<i; j++)
      {
        if(comp_hd(kmer_final[i].kmer, kmer_final[j].kmer)<=hd_thr)
        {
          psm_set[i].set.push_back(j);
          psm_set[j].set.push_back(i);
        }
      }
      
      for(int j=top_num; j<kmer_final.size(); j++)
      {
        if(comp_hd(kmer_final[i].kmer, kmer_final[j].kmer)<=hd_thr)
        {
          psm_set[i].set.push_back(j);
        }
      }
    }
    
    for(int i=0; i<psm_set.size(); i++)
    {
      psm_set[i].occurr=0;
      
      for(int j=0; j<psm_set[i].set.size(); j++)
      {
        psm_set[i].occurr+=kmer_final[psm_set[i].set[j]].occurr;
      }
    }
    
    cout<<"Finished constructing PSMs.\n";
    cout<<"There are altogether "<<psm_set.size()<<" PSMs constructed.\n";
  }
  
  
  // Compute the z-scores of all PSMs
  void test_psm(vector<psm>& psm_set, float& sum1, float& sum2)
  {
    for(int i=0; i<psm_set.size(); i++)
    {
      psm_set[i].score=ztest(psm_set[i].cover1, psm_set[i].cover2, sum1, sum2);
    }
    
    
    cout<<"Finished computing the z-scores for all PSMs\n";
  }
  
  
  // Initialize hte PSM amtrices
  void init_mat(vector<float_vector>& psm_mat, int& length)
  {
    for(int i=0; i<psm_mat.size(); i++)
    {
      psm_mat[i].clear();
    }
    psm_mat.clear();
    
    vector<float> tmp_vec;
    
    for(int i=0; i<4; i++)
    {
      tmp_vec.clear();
      for(int j=0; j<length; j++)
      {
        tmp_vec.push_back(0);
      }
      psm_mat.push_back(tmp_vec);
    }
  }
  
  
  // Compute the PSM matrix based on k-mer set of PSM
  void comp_mat(vector<float_vector>& psm_mat, vector<int>& psm_set, vector<kmer_all>& kmer_final, int& length)
  {
    int i, j;
    
    init_mat(psm_mat, length);
    
    for(i=0; i<psm_set.size(); i++)
    {
      for(j=0; j<length; j++)
      {
        psm_mat[r_alphabet[kmer_final[psm_set[i]].kmer[j]]][j]+=kmer_final[psm_set[i]].occurr;
      }
    }
  }
  
  
  // Normalize the matrices of the PSMs
  void norm_mat(vector<float_vector>& tmp_mat, int& occurr, int& length)
  {
    for(int i=0; i<4; i++)
    {
      for(int j=0; j<length; j++)
      {
        tmp_mat[i][j]/=occurr;
      }
    }
  }
  
  
  // Compute the matrices for all PSMs
  void fill_psm(vector<psm>& psm_set, vector<kmer_all>& kmer_final, int& length)//, map<char, int>& r_alphabet)
  {
    int i;
    
    for(i=0; i<psm_set.size(); i++)
    {
      comp_mat(psm_set[i].mat, psm_set[i].set, kmer_final, length);
      norm_mat(psm_set[i].mat, psm_set[i].occurr, length);
    }
    
    cout<<"Finished computing the matrices information of each PSM.\n";
  }
  
  
  // Initialize the information of the consensus string
  void init_str(string& str, int length)
  {
    str.assign(length, 'A');
  }
  
  
  // Choose the element with the maximum value in this column
  int max_col(vector<float_vector>& mat, int col)
  {
    float max_val=0;
    int max_id=0;
    
    for(int i=0; i<mat.size(); i++)
    {
      if(mat[i][col]>max_val)
      {
        max_val=mat[i][col];
        max_id=i;
      }
    }
    
    return max_id;
  }
  
  
  // Compute the consensus string for the current PSM
  void get_cons(string& cons, vector<float_vector>& mat)
  {
    int i;
    
    init_str(cons, mat[0].size());
    
    for(i=0; i<cons.length(); i++)
    {
      cons[i]=alphabet[max_col(mat, i)];
    }
  }
  
  
  // Compute the consensus strings of the PSMs
  void comp_cons(vector<psm>& psm_set)
  {
    int i;
    
    for(i=0; i<psm_set.size(); i++)
    {
      get_cons(psm_set[i].cons, psm_set[i].mat);
    }
    
    cout<<"Finished computing the consensus strings for each PSM.\n";
  }
  
  
  // Add one flanking segment to the matrix
  void add_str(vector<float_vector>& tmp_mat, string str)
  {
    for(int i=0; i<str.size(); i++)
    {
      tmp_mat[r_alphabet[str[i]]][i]++;
    }
  }
  
  
  // get the left and right flanking matrices
  void get_flank(vector<float_vector>& tmp_mat, vector<int>& site, sequence& seq_final, int direct, int id)
  {
    init_mat(tmp_mat, lmer_length);
    
    for(int i=0; i<site.size(); i+=2)
    {
      if(direct==-1)
      {
        add_str(tmp_mat, seq_final.strand[site[i]].substr(site[i+1]-lmer_length, lmer_length));
      }
      else if(direct==1)
      {
        add_str(tmp_mat, seq_final.strand[site[i]].substr(site[i+1]+kmer_length, lmer_length));
      }
      else
      {
        cerr<<"Error: There's errors in function get_flank. Please check the parameters.\n";
        exit(0);
      }
    }
  }
  
  
  // Combine two sets of k-mers
  int combine_major(kmer_set& major_set, kmer_set& minor_set)
  {
    int res=major_set.kmer.size();
    major_set.kmer.insert(major_set.kmer.end(), minor_set.kmer.begin(), minor_set.kmer.end());
    major_set.occurr.insert(major_set.occurr.end(), minor_set.occurr.begin(), minor_set.occurr.end());
    major_set.cover1.insert(major_set.cover1.end(), minor_set.cover1.begin(), minor_set.cover1.end());
    major_set.cover2.insert(major_set.cover2.end(), minor_set.cover2.begin(), minor_set.cover2.end());
    major_set.score.insert(major_set.score.end(), minor_set.score.begin(), minor_set.score.end());
    major_set.site.insert(major_set.site.end(), minor_set.site.begin(), minor_set.site.end());
    
    return res;
  }
  
  
  // Obtain the k-mer information as well as their flanking segments
  void fill_kmer(kmer_set& major_set, vector<kmer_all>& kmer_final, int seq_num, sequence& seq_final)
  {
    kmer_all tmp_kmer;
    
    for(int i=0; i<major_set.kmer.size(); i++)
    {
      tmp_kmer.kmer=major_set.kmer[i];
      tmp_kmer.occurr=major_set.occurr[i];
      tmp_kmer.cover1=major_set.cover1[i];
      tmp_kmer.cover2=major_set.cover2[i];
      tmp_kmer.site=major_set.site[i];
      tmp_kmer.score=major_set.score[i];
      get_flank(tmp_kmer.left, major_set.site[i], seq_final, -1, i);
      get_flank(tmp_kmer.right, major_set.site[i], seq_final, 1, i);
      kmer_final.push_back(tmp_kmer);
    }
    
    cout<<"Finished calculating the profile matrices by the flanking segments beside the "<<kmer_length<<"-mers.\n";
  }
  
  
  // Free the memory of kmer_set
  void free_kmer(kmer_set& major_set)
  {
    for(int i=0; i<major_set.kmer.size(); i++)
    {
      string().swap(major_set.kmer[i]);
      vector<int>().swap(major_set.site[i]);
    }
    vector<string>().swap(major_set.kmer);
    vector<int_vector>().swap(major_set.site);
  }
  
  
  // Check whether two vectors have overlap
  float comp_sw(vector<float_vector>& mat1, vector<float_vector>& mat2)
  {
    if(mat1[0].size() != mat2[0].size())
    {
      cerr<<"Error: The two matrices are not equal in width!\n";
      exit(1);
    }
    
    int length = mat1[0].size();
    float sw=0;
    
    for(int j = 0; j < length; j++)
    {
      for(int i=0; i<4; i++)
      {
        sw+=(mat1[i][j] - mat2[i][j]) * 
          (mat1[i][j] - mat2[i][j]);
      }
    }
    
    sw/=length;
    sw=2-sw;
    
    return sw;
  }
  
  
  // check whether two PSMs have overlap
  float comp_overlap(vector<int>& set1, vector<int>& set2)
  {
    int flag=0;
    
    for(int i=0; i<set1.size(); i++)
    {
      if(find(set2.begin(), set2.end(), set1[i])
           !=set2.end() && set1[i]<top_num)
      {
        flag=1;
        break;
      }
    }
    
    return flag;
  }
  
  
  // build similarity graph of PSMs
  void build_graph(vector<psm>& psm_set, vector<int_vector>& psm_nbr)
  {
    vector<int> tmp_vec;
    
    for(int i=0; i<psm_set.size(); i++)
    {
      psm_nbr.push_back(tmp_vec);
    }
    
    for(int i=1; i<psm_set.size(); i++)
    {
      for(int j=0; j<i; j++)
      {
        float tmp_sw=comp_sw(psm_set[i].mat, psm_set[j].mat);
        if(tmp_sw>=sw_thr)
        {
          psm_nbr[i].push_back(j);
          psm_nbr[j].push_back(i);
        }
      }
    }
    
    cout<<"Finished building PSM Adjacency Graph.\n";
    cout<<"There're altogether "<<psm_nbr.size()<<" nodes in this graph.\n";
  }
  
  
  // transform PSM set into k-mer set
  void psm2kmer(vector<int>& tmp_kmer, vector<int>& tmp_psm, vector<psm>& all_psm, int kmer_num)
  {
    tmp_kmer.clear();
    vector<int> flag;
    generate_flag(kmer_num, flag);
    
    for(int i=0;i<tmp_psm.size(); i++)
    {
      mask_kmer(all_psm[tmp_psm[i]].set, flag);
    }
    
    for(int i=0; i<kmer_num; i++)
    {
      if(flag[i]==1)
      {
        tmp_kmer.push_back(i);
      }
    }
  }
  
  
  // calculate the approximate coverage of the current k-mer set
  float comp_cover(vector<kmer_all>& kmer_final, vector<int>& tmp_set, int flag)
  {
    int tmp_cover=0;
    
    if(flag==1)
    {
      for(int i=0; i<tmp_set.size(); i++)
      {
        tmp_cover+=kmer_final[tmp_set[i]].cover1;
      }
    }
    else if(flag==-1)
    {
      for(int i=0; i<tmp_set.size(); i++)
      {
        tmp_cover+=kmer_final[tmp_set[i]].cover2;
      }
    }
    
    return tmp_cover;
  }
  
  
  // compute the coverage of PSMs
  void cover_psm(vector<psm>& psm_set, vector<kmer_all>& kmer_final)
  {
    for(int i=0; i<psm_set.size(); i++)
    {
      psm_set[i].cover1=comp_cover(kmer_final, psm_set[i].set, 1);
      psm_set[i].cover2=comp_cover(kmer_final, psm_set[i].set, -1);
    }
    
    cout<<"Finished cumputing the coverage for all PSMs.\n";
  }
  
  
  // calculate the occurrence of the current k-mer set
  int comp_occurr(vector<kmer_all>& kmer_final, vector<int>& tmp_set)
  {
    int tmp_occurr=0;
    
    for(int i=0; i<tmp_set.size(); i++)
    {
      tmp_occurr+=kmer_final[tmp_set[i]].occurr;
    }
    
    return tmp_occurr;
  }
  
  
  // calculate the SW similarity score between two columns
  float sw_col(vector<float_vector>& mat1, vector<float_vector>& mat2, int col_id)
  {
    float result=2;
    
    for(int i=0; i<4; i++)
    {
      result-=(mat1[i][col_id]-mat2[i][col_id])*(mat1[i][col_id]-mat2[i][col_id]);
    }
    
    return result;
  }
  
  
  // calculate the SW similarity scores
  float sw_sim(vector<float_vector>& mat1, vector<float_vector>& mat2)
  {
    float result=0;
    
    for(int j=0; j<mat1[0].size(); j++)
    {
      result+=sw_col(mat1, mat2, j);
    }
    
    return result;
  }
  
  
  // Get one random index based on one distribution
  int get_rand(vector<float>& sim_vec)
  {
    float* p;
    p=new float [sim_vec.size()];
    
    for(int i=0; i<sim_vec.size(); i++)
    {
      p[i]=sim_vec[i];
    }
    
    int result=genmulone(p, (long)sim_vec.size());
    
    free_pt(p);
    
    return result;
  }
  
  
  /*************************************************************************************************
   *
   * Function sam_node used to choose the new node to initiate the Gibbs Sampling randomly
   *
   ************************************************************************************************/
  
  int sam_node(vector<int>& tmp_pre, vector<float_vector>& tmp_mat, vector<psm>& psm_set)
  {
    // int result;
    
    vector<float> sim_vec;	// vector to hold the similarity scores
    sim_vec.reserve(tmp_pre.size());
    //	float* sim_vec; // probability vector
    //  sim_vec = new float [tmp_pre.size()]; // set size
    
    for(int i=0; i<tmp_pre.size(); i++)
    {
      sim_vec.push_back(sw_sim(tmp_mat, psm_set[tmp_pre[i]].mat));
    }
    
    int result = get_rand(sim_vec); 
    // result=0;
    // float max_sim=sim_vec[0];
    // for(int i=1; i<sim_vec.size(); i++)
    // {
    // 	if(max_sim<sim_vec[i])
    // 	{
    // 		max_sim=sim_vec[i];
    // 		result=i;
    // 	}
    // }
    
    return tmp_pre[result];
  }
  
  
  /*************************************************************************************************
   *
   * Function choose_node used to choose the old node to initiate the Gibbs Sampling uniformly
   *
   ************************************************************************************************/
  
  int choose_node(vector<int>& tmp_pre)
  {
    float* sam_vec;
    sam_vec=new float [tmp_pre.size()];
    
    for(int i=0; i<tmp_pre.size(); i++)
    {
      sam_vec[i]=1/(float)tmp_pre.size();
    }
    
    int res;
    res=genmulone(sam_vec, tmp_pre.size());
    res=tmp_pre[res];
    
    free_pt(sam_vec);
    return res;
  }
  
  
  /*************************************************************************************************
   *
   * Function get_nbr used to get the neighbor list of old node
   *
   ************************************************************************************************/
  
  void get_nbr(int& old_node, vector<int>& psm_flag, vector<int>& tmp_nbr, vector<int_vector>& psm_nbr, vector<int>& tmp_pre)
  {
    tmp_nbr.clear();
    
    for(int i=0; i<psm_nbr[old_node].size(); i++)
    {
      if(psm_flag[psm_nbr[old_node][i]]==0 && find(tmp_pre.begin(), tmp_pre.end(), psm_nbr[old_node][i])!=tmp_pre.end())
      {
        tmp_nbr.push_back(psm_nbr[old_node][i]);
      }
    }
  }
  
  
  /*************************************************************************************************
   *
   * Function comp_sc used to calculate the motif score of the current preliminary motif
   *
   ************************************************************************************************/
  
  float comp_sc(vector<float_vector> mat, int& occurr, int& length)
  {
    float sc=0;
    
    for(int i=0; i<4; i++)
    {
      for(int j=0; j<length; j++)
      {
        if(mat[i][j]==0)
        {
          mat[i][j]=0.00001;
        }
        
        mat[i][j]=mat[i][j]*log(mat[i][j]/nt_freq[i]);
        sc+=mat[i][j];
      }
    }
    
    sc/=length;
    sc=exp(sc);
    sc=sc*occurr;
    
    return sc;
  }
  
  
  /*************************************************************************************************
   *
   * Function del_element used to delete element in vector
   *
   ************************************************************************************************/
  
  void del_element(vector<int>& tmp_vec, int& val)
  {
    for(vector<int>::iterator i=tmp_vec.begin(); i!=tmp_vec.end(); i++)
    {
      if(*i==val)
      {
        tmp_vec.erase(i);
        break;
      }
    }
  }
  
  
  /*************************************************************************************************
   *
   * Function if_accept used to check whether we add one node to the current preliminary motif
   *
   ************************************************************************************************/
  
  int if_accept(vector<int>& tmp_pre, int& tmp_occurr, vector<float_vector>& tmp_mat, int& old_node, 
                int& new_node, vector<kmer_all>& kmer_final, vector<psm>& psm_set, int seq_num)
  {
    vector<int> old_pre=tmp_pre;
    old_pre.push_back(old_node);
    vector<int> old_kmer;
    psm2kmer(old_kmer, old_pre, psm_set, kmer_final.size());
    vector<float_vector> old_mat;
    int old_occurr=comp_occurr(kmer_final, old_kmer);
    comp_mat(old_mat, old_kmer, kmer_final, kmer_length);
    norm_mat(old_mat, old_occurr, kmer_length);
    float old_score=comp_sc(old_mat, old_occurr, kmer_length);
    
    vector<int> new_pre=old_pre;
    new_pre.push_back(old_node);
    new_pre.push_back(new_node);
    vector<int> new_kmer;
    psm2kmer(new_kmer, new_pre, psm_set, kmer_final.size());
    int new_occurr=comp_occurr(kmer_final, new_kmer);
    vector<float_vector> new_mat;
    comp_mat(new_mat, new_kmer, kmer_final, kmer_length);
    norm_mat(new_mat, new_occurr, kmer_length);
    float new_score=comp_sc(new_mat, new_occurr, kmer_length);
    
    vector<int> del_pre=old_pre;
    del_pre.push_back(new_node);
    vector<int> del_kmer;
    psm2kmer(del_kmer, del_pre, psm_set, kmer_final.size());
    int del_occurr=comp_occurr(kmer_final, del_kmer);
    vector<float_vector> del_mat;
    comp_mat(del_mat, del_kmer, kmer_final, kmer_length);
    norm_mat(del_mat, del_occurr, kmer_length);
    float del_score=comp_sc(del_mat, del_occurr, kmer_length);
    
    vector<float> v_score; 
    v_score.push_back(old_score);
    v_score.push_back(del_score);
    v_score.push_back(new_score);
    
    double sum = accumulate(v_score.begin(), v_score.end(), 0); // get the sum of the vector
    
    for (int i = 0; i < v_score.size(); i ++)
    {
      v_score[i] /= sum;
    }
    
    return(get_rand(v_score)); // randomly sample an id according to the probability distributation
    // calculated based on the motif score
  }
  
  
  /*************************************************************************************************
   *
   * Function psm2pre used to transform PSMs into preliminary motifs
   *
   ************************************************************************************************/
  
  void psm2pre(vector<int_vector>& pre_out, vector<psm>& psm_set)
  {
    pre_out.clear();
    
    for(int i=0; i<psm_set.size(); i++)
    {
      pre_out.push_back(psm_set[i].set);
    }
  }
  
  
  /*************************************************************************************************
   *
   * Function gibbs_sam used to perform Gibbs Sampling
   *
   ************************************************************************************************/
  
  void gibbs_sam(vector<int_vector>& pre_out, vector<psm>& psm_set, 
                 vector<kmer_all>& kmer_final, int& num_pre, int& num_iter, 
                 vector<int_vector>& psm_nbr, int seq_num)
  {
    vector<int> psm_flag;
    generate_flag(psm_set.size(), psm_flag);
    int seed_psm;
    vector<int> tmp_pre, tmp_kmer, tmp_nbr;
    vector<float_vector> tmp_mat;
    int tmp_occurr, old_node, new_node, del_node;
    float tmp_score, tmp_cover1, tmp_cover2;
    int sam_flag;
    
    cout<<"Begin Gibbs Sampling ...\n";
    
    for(int i=1; i<=num_pre; i++)
    {
      seed_psm=find_seed(psm_flag);
      
      if(seed_psm==-1)
      {
        break;
      }
      
      tmp_pre.clear();
      tmp_pre.push_back(seed_psm);
      for(int j=0; j<psm_nbr[seed_psm].size(); j++)
      {
        if(psm_flag[psm_nbr[seed_psm][j]]==0)
        {
          tmp_pre.push_back(psm_nbr[seed_psm][j]);
        }
      }
      mask_kmer(tmp_pre, psm_flag);
      psm2kmer(tmp_kmer, tmp_pre, psm_set, kmer_final.size());
      
      if(tmp_pre.size()==1)
      {
        pre_out.push_back(tmp_kmer);
        continue;
      }
      
      for(int j=1; j<=num_iter; j++)
      {
        old_node=choose_node(tmp_pre);
        del_element(tmp_pre, old_node);
        psm_flag[old_node]=0;
        psm2kmer(tmp_kmer, tmp_pre, psm_set, kmer_final.size());
        tmp_occurr=comp_occurr(kmer_final, tmp_kmer);
        comp_mat(tmp_mat, tmp_kmer, kmer_final, kmer_length);
        norm_mat(tmp_mat, tmp_occurr, kmer_length);
        tmp_score=comp_sc(tmp_mat, tmp_occurr, kmer_length);
        get_nbr(old_node, psm_flag, tmp_nbr, psm_nbr, tmp_pre);
        
        if(tmp_nbr.size() == 0)
        {
          psm_flag[old_node]=1;
          tmp_pre.push_back(old_node);
          continue;
        }
        
        new_node = sam_node(tmp_nbr, tmp_mat, psm_set); // select node randomly
        sam_flag = if_accept(tmp_pre, tmp_occurr, tmp_mat, old_node, 
                             new_node, kmer_final, psm_set, seq_num);
        
        if(sam_flag == 2)
        {
          tmp_pre.push_back(old_node);
          tmp_pre.push_back(new_node);
        }
        else if(sam_flag == 1)
        {
          tmp_pre.push_back(new_node);
        }
        else
        {
          tmp_pre.push_back(old_node);
        }
        
        mask_kmer(tmp_pre, psm_flag);
      }
      
      pre_out.push_back(tmp_kmer);
    }
    
    cout<<"Finished Gibbs Sampling.\n" << "We have obtained altogether " << pre_out.size() << 
      " preliminary motifs."<<endl;
  }
  
  
  /*************************************************************************************************
   *
   * Function build_flank used to build flanking matrices
   *
   ************************************************************************************************/
  
  void build_flank(vector<int>& tmp_set, vector<kmer_all>& kmer_final, vector<float_vector>& flank, int direct, int length)
  {
    init_mat(flank, length);
    
    for(int i=0; i<tmp_set.size(); i++)
    {
      if(direct==-1)
      {
        for(int j=0; j<4; j++)
        {
          for(int k=0; k<length; k++)
          {
            flank[j][k]+=kmer_final[tmp_set[i]].left[j][k];
          }
        }
      }
      else
      {
        for(int j=0; j<4; j++)
        {
          for(int k=0; k<length; k++)
          {
            flank[j][k]+=kmer_final[tmp_set[i]].right[j][k];
          }
        }
      }
    }
  }
  
  
  /*************************************************************************************************
   *
   * Function site_pre used to combine the site information of all k-mers contained
   *
   ************************************************************************************************/
  
  void site_pre(pre_mtf& tmp_pre, vector<kmer_all>& kmer_final)
  {
    tmp_pre.site.clear();
    pair<int, int> tmp_pair;
    map<pair<int, int>, int> site_map;
    
    for(int i=0; i<tmp_pre.set.size(); i++)
    {
      for(int j=0; j<kmer_final[tmp_pre.set[i]].site.size(); j+=2)
      {
        tmp_pair=make_pair(kmer_final[tmp_pre.set[i]].site[j], 
                           kmer_final[tmp_pre.set[i]].site[j+1]);
        site_map[tmp_pair]=1;
      }
    }
    
    for(map<pair<int, int>, int>::iterator it=site_map.begin(); it!=site_map.end(); it++)
    {
      tmp_pre.site.push_back((it->first).first);
      tmp_pre.site.push_back((it->first).second);
    }
  }
  
  
  /*************************************************************************************************
   *
   * Function get_pre used to get the set of preliminary motifs based on preliminary output 
   * from Gibbs Sampling
   *
   ************************************************************************************************/
  
  void get_pre(vector<pre_mtf>& all_pre, vector<int_vector>& pre_out, 
               vector<kmer_all>& kmer_final, int seq_num)
  {
    for(int i=0; i<pre_out.size(); i++)
    {
      pre_mtf tmp_pre;
      
      tmp_pre.set=pre_out[i];
      site_pre(tmp_pre, kmer_final);
      tmp_pre.occurr=comp_occurr(kmer_final, tmp_pre.set);
      comp_mat(tmp_pre.mat, tmp_pre.set, kmer_final, kmer_length);
      norm_mat(tmp_pre.mat, tmp_pre.occurr, kmer_length);
      
      // cout << tmp_pre.occurr << "\n";
      // exit(0);
      // for (int j = 0; j < tmp_pre.mat.size(); j ++)
      // {
      //   for (int k = 0; k < tmp_pre.mat[j].size(); k ++)
      //   {
      //     cout << tmp_pre.mat[j][k] << " ";
      //   }
      //   
      //   cout << "\n";
      // }
      
      // exit(0);
      
      build_flank(tmp_pre.set, kmer_final, tmp_pre.left, -1, lmer_length);
      build_flank(tmp_pre.set, kmer_final, tmp_pre.right, 1, lmer_length);
      norm_mat(tmp_pre.left, tmp_pre.occurr, lmer_length);
      norm_mat(tmp_pre.right, tmp_pre.occurr, lmer_length);
      tmp_pre.cover1=comp_cover(kmer_final, tmp_pre.set, 1);
      tmp_pre.cover2=comp_cover(kmer_final, tmp_pre.set, -1);
      tmp_pre.score=ztest(tmp_pre.cover1, tmp_pre.cover2, sum_exp, sum_bg);
      
      // cout << tmp_pre.cover1 << " " << tmp_pre.cover2 << " " << tmp_pre.score << "\n";
      // exit(0);
      
      all_pre.push_back(tmp_pre);
    }
  }
  
  
  /*************************************************************************************************
   *
   * Function test_cell used to check one cell in one column (Two Proportion z-test)
   *
   ************************************************************************************************/
  
  float test_cell(float this_val, float that_val, int occurr)
  {
    float pos_occurr=occurr*this_val;
    float neg_occurr=occurr*that_val;
    float p=(pos_occurr+neg_occurr)/(2*occurr);
    float z=fabs((this_val-that_val)/sqrt(p*(1-p)*(2/(float)occurr)));
    
    return z;
  }
  
  
  /*************************************************************************************************
   *
   * Function test_col used to check one column (Two Proportion z-test)
   *
   ************************************************************************************************/
  
  float test_col(vector<float_vector>& flank, int col_id, int occurr)
  {
    float max_z=0;
    float tmp_z;
    
    for(int i=0; i<4; i++)
    {
      tmp_z=test_cell(flank[i][col_id], nt_freq[i], occurr);
      
      if(tmp_z>max_z)
      {
        max_z=tmp_z;
      }
    }
    
    return max_z;
  }
  
  
  /*************************************************************************************************
   *
   * Function refine_one used to refine one preliminary motif
   *
   ************************************************************************************************/
  
  void refine_one(pre_mtf& one_pre)
  {
    if(lmer_length <= 0)
    {
      one_pre.begin = 0;
      
      return;
    }
    
    one_pre.begin=lmer_length-1;
    
    for(int i=lmer_length-1; i>=0; i--)
    {
      if(test_col(one_pre.left, i, one_pre.occurr)<thr2)
      {
        one_pre.begin=i+1;
        break;
      }
    }
    
    one_pre.end=0;
    
    for(int i=0; i<lmer_length; i++)
    {
      if(test_col(one_pre.right, i, one_pre.occurr)<thr2)
      {
        one_pre.end=i-1;
        break;
      }
    }
    
    // cout << one_pre.begin << " " << one_pre.end << "\n";
  }
  
  
  /*************************************************************************************************
   *
   * Function refine_mtf used to refine the preliminary motifs
   *
   ************************************************************************************************/
  
  void refine_mtf(vector<pre_mtf>& all_pre)
  {
    for(int i=0; i<all_pre.size(); i++)
    {
      refine_one(all_pre[i]);
      
      //		cout << all_pre[i].begin << " " << 
      //			all_pre[i].end << "\n";
    }
  }
  
  
  // Used in compute_chi_square
  static double igf(double S, double Z)
  {
    if(Z < 0.0)
    {
      return 0.0;
    }
    double Sc = (1.0 / S);
    Sc *= pow(Z, S);
    Sc *= exp(-Z);
    
    double Sum = 1.0;
    double Nom = 1.0;
    double Denom = 1.0;
    
    for(int I = 0; I < 200; I++)
    {
      Nom *= Z;
      S++;
      Denom *= S;
      Sum += (Nom / Denom);
    }
    
    return Sum * Sc;
  }
  
  
  // Used in compute_chi_square
  double gamma2(double N)
  {
    const long double SQRT2PI = 2.5066282746310005024157652848110452530069867406099383;
    
    long double Z = (long double)N;
    long double Sc = powl((Z + A), (Z + 0.5));
    Sc *= expl(-1.0 * (Z + A));
    Sc /= Z;
    
    long double F = 1.0;
    long double Ck;
    long double Sum = SQRT2PI;
    
    
    for(int K = 1; K < A; K++)
    {
      Z++;
      Ck = powl(A - K, K - 0.5);
      Ck *= expl(A - K);
      Ck /= F;
      
      Sum += (Ck / Z);
      
      F *= (-1.0 * K);
    }
    
    return (double)(Sum * Sc);
  }
  
  
  /*************************************************************************************************
   *
   * Function compute_chi_square used to calculate the chi square p-value
   * Link: https://www.codeproject.com/Articles/432194/How-to-Calculate-the-Chi-Squared-P-Value
   *
   ************************************************************************************************/
  
  double compute_chi_square(int Dof, double Cv)
  {
    if(Cv < 0 || Dof < 1)
    {
      return 0.0;
    }
    double K = ((double)Dof) * 0.5;
    double X = Cv * 0.5;
    if(Dof == 2)
    {
      return exp(-1.0 * X);
    }
    
    double PValue = igf(K, X);
    if(isnan(PValue) || isinf(PValue) || PValue <= 1e-8)
    {
      return 1e-14;
    } 
    
    PValue /= gamma2(K);
    //PValue /= tgamma(K); 
    
    return (1.0 - PValue);
  }
  
  
  /*************************************************************************************************
   *
   * Function compute_pvalue used to sift preliminary motifs using chi square p-values
   *
   ************************************************************************************************/
  
  void compute_pvalue(vector<pre_mtf>& all_pre, vector<kmer_all>& kmer_final, 
                      map<string, float>& kmer_bg)
  {
    vector<pre_mtf> out_pre; // declare a vector to hold retained preliminary motifs
    out_pre.reserve(all_pre.size()); // reserve space

    for (int i = 0; i < all_pre.size(); i ++) // for each preliminary motif
    {
      //cout << i << "\n";
      // int n = 0; // initialize the counting number
      // calculate the chi square p-value
      double chi_val = 0; // initialize the chi-square static
      int kmer_count = 0;
      float other_kmer = 1; 
      
      // cout << "Size: " << all_pre[i].set.size() << "\n";
      
      for (int j = 0; j < all_pre[i].set.size(); j ++) // for each k-mer in the preliminary motif
      {
        // n ++; // update the number of k-mers
        // double P = T.AA[kmer_final[all_pre[i].set[j]].kmer.substr(0, 1)] * 
        //   T.B[kmer_final[all_pre[i].set[j]].kmer.substr(0, 2)] * 
        //   T.C[kmer_final[all_pre[i].set[j]].kmer.substr(0, 3)]; // probability of the first three elements
        // // in the string
        // 
        // // cout << "probability: " << P << "\n";
        // // exit(0);
        // 
        // for (int k = 3; k < kmer_final[all_pre[i].set[j]].kmer.length(); k ++) // for the remaining elements
        // {
        //   P *= T.D[kmer_final[all_pre[i].set[j]].kmer.substr(k - 3, 4)]; // probabilities of 
        //   // the remaining elements
        //   // maybe I can logrithmize the p-value for acceleration
        // }
        double P = kmer_bg[kmer_final[all_pre[i].set[j]].kmer]; //k-mer frequency in background
        chi_val += pow((kmer_final[all_pre[i].set[j]].occurr / all_pre[i].occurr - P), 2) / P;
        // calculate the chi-square statistic
        other_kmer -= P;
        kmer_count += kmer_final[all_pre[i].set[j]].occurr; 
      }
      chi_val += pow(((all_pre[i].occurr - kmer_count) / all_pre[i].occurr - other_kmer), 2) / 
        other_kmer;
      
      // cout << chi_val << " " << n << "\n";
      all_pre[i].pvalue = compute_chi_square(all_pre[i].set.size(), chi_val); // calculate the p-value based on
      // the chi square value
      
      // cout << chi_val << " " << all_pre[i].pvalue << "\n";
    }
    
    // p-value adjustment
    for (int i = 0; i < all_pre.size(); i ++)
    {
      all_pre[i].padj = all_pre[i].pvalue * all_pre.size();
      
      // cout << "p-adj: " << all_pre[i].padj << "\n";
    }
  }
  
  
  /*************************************************************************************************
   *
   * Function cut_mat used to transform matrices
   *
   ************************************************************************************************/
  
  void cut_mat(pre_mtf& tmp_pre, vector<float_vector>& mat)
  {
    mat.clear();
    deque<float> tmp_ln;
    vector<float> tmp_vec;
    
    for (int i = 0; i < 4; i ++)
    {
      tmp_ln.clear();
      tmp_vec.clear();
      
      for (int j = 0; j < kmer_length; j ++)
      {
        tmp_ln.push_back
        (tmp_pre.mat[i][j]);
      }
      
      if (lmer_length > 0)
      {
        if (tmp_pre.end >= 0)
        {
          for (int j = 0; j <= tmp_pre.end; j ++)
          {
            tmp_ln.push_back
            (tmp_pre.right[i][j]);
          }
        }
        
        if (tmp_pre.begin < lmer_length)
        {
          for (int j = lmer_length - 1; 
               j >= tmp_pre.begin; j --)
          {
            tmp_ln.push_front
            (tmp_pre.left[i][j]);
          }
        }
      }
      
      for (deque<float>::iterator i = tmp_ln.begin(); 
           i != tmp_ln.end(); i ++)
      {
        tmp_vec.push_back(*i);
      }
      
      mat.push_back(tmp_vec);
    }
  }
  
  
  /*************************************************************************************************
   *
   * Function make_wild used to make hash table of wild cards
   *
   ************************************************************************************************/
  
  void make_wild(char* wild_card, char* reverse_wild)
  {
    wild_card[0]='N';
    reverse_wild[0]='N';
    wild_card[1]='A';
    reverse_wild[1]='T';
    wild_card[2]='C';
    reverse_wild[2]='G';
    wild_card[3]='M';	// A, C
    reverse_wild[3]='K';
    wild_card[4]='G';
    reverse_wild[4]='C';
    wild_card[5]='R';	// A, G
    reverse_wild[5]='Y';
    wild_card[6]='S';	// C, G
    reverse_wild[6]='S';
    wild_card[7]='V';	// A, C, G
    reverse_wild[7]='B';
    wild_card[8]='T';
    reverse_wild[8]='A';
    wild_card[9]='W';	// A, T
    reverse_wild[9]='W';
    wild_card[10]='Y';	// C, T
    reverse_wild[10]='R';
    wild_card[11]='H';	// A, C, T
    reverse_wild[11]='D';
    wild_card[12]='K';	// G, T
    reverse_wild[12]='M';
    wild_card[13]='D';	// A, G, T
    reverse_wild[13]='H';
    wild_card[14]='B';	// G, C, T
    reverse_wild[14]='V';
    
    cout<<"Finished constructing wild card tables.\n";
  }
  
  
  /*************************************************************************************************
   *
   * Function col_deg used to get the degenerate consensus string 
   *
   ************************************************************************************************/
  
  void col_deg(string& deg, string& rev_deg, vector<float_vector>& mat, int& col_id)
  {
    int sum=0;
    
    for(int i=0; i<4; i++)
    {
      if(mat[i][col_id]>0.25)
      {
        sum+=(int)(rint(pow(2, i)));
      }
    }
    
    deg[col_id]=wild_card[sum];
    rev_deg[rev_deg.length()-1-col_id]=reverse_wild[sum];
  }
  
  
  /*************************************************************************************************
   *
   * Function get_deg used to get the degenerate consensus strings
   *
   ************************************************************************************************/
  
  void get_deg(string& deg, string& rev_deg, vector<float_vector>& mat)
  {
    init_str(deg, mat[0].size());
    init_str(rev_deg, mat[0].size());
    
    for(int i=0; i<deg.length(); i++)
    {
      col_deg(deg, rev_deg, mat, i);
    }
  }
  
  
  /*************************************************************************************************
   *
   * Function build_rev used to get the reverse consensus string
   *
   ************************************************************************************************/
  
  void build_rev(string& rev_cons, string& cons)
  {
    init_str(rev_cons, cons.length());
    
    for(int i=0; i<cons.length(); i++)
    {
      rev_cons[cons.length()-1-i]=alphabet[3-r_alphabet[cons[i]]];
    }
  }
  
  
  /*************************************************************************************************
   *
   * Function build_mtf used to transform one preliminary motif into one motif
   *
   ************************************************************************************************/
  
  void build_mtf(pre_mtf& tmp_pre, mtf& tmp_mtf)
  {
    tmp_mtf.set=tmp_pre.set;
    tmp_mtf.nsites=tmp_pre.occurr;
    cut_mat(tmp_pre, tmp_mtf.mat);
    tmp_mtf.alength=tmp_mtf.mat[0].size();
    get_cons(tmp_mtf.cons, tmp_mtf.mat);
    get_deg(tmp_mtf.deg, tmp_mtf.rev_deg, tmp_mtf.mat);
    build_rev(tmp_mtf.rev_cons, tmp_mtf.cons);
    tmp_mtf.score=tmp_pre.score;
  }
  
  
  /*************************************************************************************************
   *
   * Function get_site used to get the sites of motif
   *
   ************************************************************************************************/
  
  void get_site(pre_mtf& tmp_pre, mtf& tmp_mtf, sequence& seq_final)
  {
    st tmp_st;
    
    if(str_flag==1)
    {
      for(int i=0; i<tmp_pre.site.size(); i+=2)
      {
        tmp_st.header=seq_final.name[tmp_pre.site[i]];
        tmp_st.strand='+';
        tmp_st.begin=tmp_pre.site[i+1]-lmer_length+tmp_pre.begin+1;
        tmp_st.end=tmp_pre.site[i+1]+kmer_length+tmp_pre.end+1;
        
        tmp_mtf.site.push_back(tmp_st);
      }
    }
    else
    {
      for(int i=0; i<tmp_pre.site.size(); i+=2)
      {
        tmp_st.header=seq_final.name[tmp_pre.site[i]/2];
        
        if(tmp_pre.site[i]%2==0)
        {
          tmp_st.strand='+';
          //				cout << tmp_pre.site[i+1] << "\n";
          tmp_st.begin=tmp_pre.site[i+1]-lmer_length+tmp_pre.begin+1;
          tmp_st.end=tmp_pre.site[i+1]+kmer_length+tmp_pre.end+1;
          //				cout << tmp_st.begin << " " << lmer_length 
          //					<< " " << tmp_pre.begin << "\n";
        }
        else
        {
          tmp_st.strand='-';
          tmp_st.begin=tmp_pre.site[i+1]-lmer_length+tmp_pre.begin;
          tmp_st.end=tmp_pre.site[i+1]+kmer_length-1+tmp_pre.end;
          int tmp_val=tmp_st.begin;
          tmp_st.begin=seq_final.strand[tmp_pre.site[i]].
          size()-tmp_st.end-1;
          tmp_st.end=seq_final.strand[tmp_pre.site[i]].size()-tmp_val;
        }
        
        tmp_mtf.site.push_back(tmp_st);
      }
    }
  }
  
  
  /*************************************************************************************************
   *
   * Function pre2mtf used to transform preliminary motifs into motifs
   *
   ************************************************************************************************/
  
  void pre2mtf(vector<pre_mtf>& all_pre, vector<mtf>& all_mtf, sequence& seq_final)
  {
    for(int i=0; i<all_pre.size(); i++)
    {
      mtf tmp_mtf;
      build_mtf(all_pre[i], tmp_mtf);
      get_site(all_pre[i], tmp_mtf, seq_final);
      all_mtf.push_back(tmp_mtf);
    }
    
    cout<<"Finished transforming preliminary motifs into motifs.\n";
  }
  
  
  /*************************************************************************************************
   *
   * Function sort_mtf used to sort all the motifs
   *
   ************************************************************************************************/
  
  void sort_mtf(vector<mtf>& mtf_set)
  {
    float *mtf_scores;
    int *mtf_id1;
    int num=mtf_set.size();
    mtf_scores=new float [num];
    mtf_id1=new int [num];
    
    int i;
    
    for(i=0; i<num; i++)
    {
      mtf_scores[i]=mtf_set[i].score;
      mtf_id1[i]=i;
    }
    
    quickSort(mtf_scores, 0, num-1, mtf_id1);
    reverse_sort(mtf_id1, num);
    
    vector<mtf> tmp_mtf;
    for(i=0; i<num; i++)
    {
      if(!last_trivial(mtf_set[mtf_id1[i]].cons))
      {
        continue;
      }
      
      tmp_mtf.push_back(mtf_set[mtf_id1[i]]);
    }
    
    mtf_set=tmp_mtf;
    
    delete [] mtf_id1;
    delete [] mtf_scores;
    
    cout << "Finished sorting the motifs according to motif scores.\n";
    cout << "Finished kicking out the trivial motifs.\nThere're " 
         << mtf_set.size() << " motifs left.\n";
  }
  
  
  /*************************************************************************************************
   *
   * Function str2vec used to transform strings into vectors
   *
   ************************************************************************************************/
  
  void str2vec(string& str, vector<char_vector>& vec)
  {
    vec.clear();
    vector<char> single_vec;
    
    for(int i=0; i<str.length(); i++)
    {
      single_vec.clear();
      
      if(str[i]=='A')
      {
        single_vec.push_back('A');
      }
      else if(str[i]=='C')
      {
        single_vec.push_back('C');
      }
      else if(str[i]=='M')
      {
        single_vec.push_back('A');
        single_vec.push_back('C');
      }
      else if(str[i]=='G')
      {
        single_vec.push_back('G');
      }
      else if(str[i]=='R')
      {
        single_vec.push_back('A');
        single_vec.push_back('G');
      }
      else if(str[i]=='S')
      {
        single_vec.push_back('C');
        single_vec.push_back('G');
      }
      else if(str[i]=='V')
      {
        single_vec.push_back('A');
        single_vec.push_back('C');
        single_vec.push_back('G');
      }
      else if(str[i]=='T')
      {
        single_vec.push_back('T');
      }
      else if(str[i]=='W')
      {
        single_vec.push_back('A');
        single_vec.push_back('T');
      }
      else if(str[i]=='Y')
      {
        single_vec.push_back('C');
        single_vec.push_back('T');
      }
      else if(str[i]=='H')
      {
        single_vec.push_back('A');
        single_vec.push_back('C');
        single_vec.push_back('T');
      }
      else if(str[i]=='K')
      {
        single_vec.push_back('G');
        single_vec.push_back('T');
      }
      else if(str[i]=='D')
      {
        single_vec.push_back('A');
        single_vec.push_back('G');
        single_vec.push_back('T');
      }
      else if(str[i]=='B')
      {
        single_vec.push_back('G');
        single_vec.push_back('C');
        single_vec.push_back('T');
      }
      else
      {
        single_vec.push_back('A');
        single_vec.push_back('C');
        single_vec.push_back('G');
        single_vec.push_back('T');
      }
      
      vec.push_back(single_vec);
    }
  }
  
  
  /*************************************************************************************************
   *
   * Function eval_hd used to calculate similarity between motifs of the same length
   *
   ************************************************************************************************/
  
  int eval_hd(string str1, string str2)
  {
    if(str1.length()!=str2.length())
    {
      cerr<<"Error: The lengths of the two strings are not equal. Please check the details.\n";
      exit(1);
    }
    
    vector<char_vector> str_vec1, str_vec2;
    str2vec(str1, str_vec1);
    str2vec(str2, str_vec2);
    
    int dist=0;
    
    for(int i=0; i<str_vec1.size(); i++)
    {
      int dt=1;
      
      for(int j=0; j<str_vec1[i].size(); j++)
      {
        if(find(str_vec2[i].begin(), str_vec2[i].end(), str_vec1[i][j])!=str_vec2[i].end())
        {
          dt=0;
          break;
        }
      }
      
      dist+=dt;
    }
    
    return dist;
  }
  
  
  /*************************************************************************************************
   *
   * Function cmp_mtf used to determine the similarity between two motifs in slack conditions
   *
   ************************************************************************************************/
  
  int cmp_mtf(mtf& mtf1, mtf& mtf2)
  {
    int length;
    if(mtf1.alength >= mtf2.alength)
    {
      length=mtf2.alength;
    }
    else
    {
      length=mtf1.alength;
    }
    int min_val=length;
    int flag;
    int tmp_val;
    
    for(int i=0; i<=mtf1.alength-length; i++)
    {
      for(int j=0; j<=mtf2.alength-length; j++)
      {
        tmp_val=eval_hd(mtf1.deg.substr(i, length), mtf2.deg.substr(j, length));
        
        if(tmp_val<min_val)
        {
          min_val=tmp_val;
        }
        
        tmp_val=eval_hd(mtf1.deg.substr(i, length), mtf2.rev_deg.substr(j, length));
        
        if(tmp_val<min_val)
        {
          min_val=tmp_val;
        }
      }
    }
    
    if(min_val>=redundant_thr)
    {
      flag=0;
    }
    else
    {
      flag=1;
    }
    
    return flag;
  }
  
  
  /*************************************************************************************************
   *
   * Function sift_mtf used to delete the similar motifs with lower motif scores
   *
   ************************************************************************************************/
  
  void sift_mtf(vector<mtf>& all_mtf)
  {
    vector<int> flag;
    generate_flag(all_mtf.size(), flag);
    
    if(all_mtf.size() <= 0)
    {
      cerr << "0 motifs were identified!\n";
      exit(0);
    }
    
    for(int i=0; i<all_mtf.size()-1; i++)
    {
      if (all_mtf[i].padj > padj_cutoff)
      {
        flag[i] = 1;
        continue;
      }
      
      for(int j=i+1; j<all_mtf.size(); j++)
      {
        if (all_mtf[j].padj > padj_cutoff)
        {
          flag[j] = 1;
          continue;
        }
        
        if(flag[j]==0 && cmp_mtf(all_mtf[i], all_mtf[j])==1)
        {
          flag[j]=1;
        }
      }
    }
    
    vector<mtf> tmp_all;
    
    for(int i=0; i<all_mtf.size(); i++)
    {
      if(flag[i]==0)
      {
        tmp_all.push_back(all_mtf[i]);
      }
    }
    
    all_mtf=tmp_all;
    
    cout<<"Finished deleting similar motifs with lower motif scores.\n"<<"There are altogether "<<
      all_mtf.size() <<" motifs left."<<endl;
  }
  
  
  /*************************************************************************************************
   *
   * Function print_mtf used to print the result into one file
   *
   ************************************************************************************************/
  
  void print_mtf(char *in_file, string& bg_file, string& out_file, vector<mtf>& all_mtf)
  {
    string tmp_out=out_file;
    tmp_out+=".meme";
    
    ofstream f_out_op(tmp_out.c_str());
    
    if(!f_out_op)
    {
      cerr<<"Error: Can't open file "<<tmp_out.c_str()<<" for output!\n";
      exit(1);
    }
    
    f_out_op<<"# MEME 4.0.0\n"<<"# Command: ./ProSampler -i "<<in_file
            <<" -b "<<bg_file.c_str()<<" -o "<<out_file.c_str()<<" -d "<<num_deg<<" -m "<<num_mtf
            <<" -f "<<num_iter<<" -k "<<kmer_length<<" -l "<<lmer_length<<" -r "
            <<redundant_thr<<" -p "<<str_flag<<" -t "<<thr1<< " -w " << thr3 << " -c "<<hd_thr<<" -z "
            <<thr2<<" -h "<<help_flag<< " -P " << padj_cutoff << "\n";
    f_out_op<<endl<<"# Begin: "<<begin_pkg<<endl;
    f_out_op<<"#   End: "<<end_pkg<<endl<<endl;
    f_out_op<<"MEME Version 4\n"<<endl;
    f_out_op<<"ALPHABET= ACGT\n"<<endl;
    
    if(str_flag==1)
    {
      f_out_op<<"Strands: +\n"<<endl;
    }
    else
    {
      f_out_op<<"Strands: + -\n"<<endl;
    }
    
    f_out_op<<"Background letter frequencies (from dataset):\n";
    f_out_op.setf(ios::fixed);
    
    for(int i=0; i<3; i++)
    {
      f_out_op<<alphabet[i]<<" "<<fixed<<setprecision(PRECISE)<<nt_freq[i]<<" ";
    }
    f_out_op<<alphabet[3]<<" "<<fixed<<setprecision(PRECISE)<<nt_freq[3]<<endl<<endl;
    
    if(num_mtf>all_mtf.size())
    {
      num_mtf=all_mtf.size();
    }
    
    for(int i=0; i<num_mtf; i++)
    {
      f_out_op<<"MOTIF "<<all_mtf[i].deg<<" "
              <<all_mtf[i].rev_deg<<" ProSampler\n"<<endl;
      f_out_op<<"letter-probability matrix: alength= 4 w= "
              <<all_mtf[i].alength<<" nsites= "<<all_mtf[i].nsites << 
      " score= " << all_mtf[i].score<< " pvalue= " << setprecision(PRECISE) << all_mtf[i].pvalue << 
        " padj= " << setprecision(PRECISE) << all_mtf[i].padj << "\n";
      
      for(int j=0; j<all_mtf[i].alength; j++)
      {
        f_out_op.setf(ios::fixed);
        for(int k=0; k<3; k++)
        {
          f_out_op<<fixed<<setprecision(PRECISE)<<all_mtf[i].mat[k][j]<<" ";
        }
        f_out_op<<fixed<<setprecision(PRECISE)<<all_mtf[i].mat[3][j]<<endl;
      }
      
      f_out_op<<endl<<endl;
    }
    
    f_out_op<<endl;
    f_out_op.unsetf(ios::fixed);
    f_out_op<<"Time "<<difftime(end_t, begin_t)<<" secs.\n";
    f_out_op.close();
    
    cout<<"Finished generating the motif file "<<tmp_out<<".\n";
    cout<<"There're altogether "<<num_mtf<<" motifs output.\n";
  }
  
  
  /*************************************************************************************************
   *
   * Function print_site used to print the result into one file
   *
   ************************************************************************************************/
  
  void print_site(char *in_file, string& bg_file, string& out_file, vector<mtf>& all_mtf)
  {
    string tmp_out=out_file;
    tmp_out+=".site";
    ofstream f_out_op(tmp_out.c_str());
    
    if(!f_out_op)
    {
      cerr<<"Error: Can't open file "<<tmp_out.c_str()<<" for output!\n";
      exit(1);
    }
    
    f_out_op<<"# MEME 4.0.0\n"<<"# Command: ./ProSampler -i "<<in_file
            <<" -b "<<bg_file.c_str()<<" -o "<<out_file.c_str()<<" -d "<<num_deg<<" -m "<<num_mtf
            <<" -f "<<num_iter<<" -k "<<kmer_length<<" -l "<<lmer_length<<" -r "
            <<redundant_thr<<" -p "<<str_flag<<" -t "<<thr1<< " -w " << thr3 << " -c "<<hd_thr<<" -z "
            <<thr2<<" -h "<<help_flag<< " -P " << padj_cutoff << "\n";
    f_out_op<<endl<<"# Begin: "<<begin_pkg<<endl;
    f_out_op<<"#   End: "<<end_pkg<<endl<<endl;
    f_out_op<<"ProSampler Version 1.0.0\n"<<endl;
    f_out_op<<"ALPHABET= ACGT\n"<<endl;
    
    if(str_flag==1)
    {
      f_out_op<<"Strands: +\n"<<endl;
    }
    else
    {
      f_out_op<<"Strands: + -\n"<<endl;
    }
    
    f_out_op<<"Background letter frequencies (from dataset):\n";
    
    if(num_mtf>all_mtf.size())
    {
      num_mtf=all_mtf.size();
    }
    
    for(int i=0; i<num_mtf; i++)
    {
      f_out_op<<"MOTIF "<<all_mtf[i].deg<<" "<<all_mtf[i].rev_deg<<" ProSampler\n"<<endl;
      f_out_op<<"letter-probability matrix: alength= 4 w= "
              <<all_mtf[i].alength<<" nsites= "<<all_mtf[i].nsites << 
      " score= " << all_mtf[i].score<< " pvalue= " << all_mtf[i].pvalue << 
        "padj= " << all_mtf[i].padj << "\n";
      for(int j=0; j<all_mtf[i].site.size(); j++)
      {
        f_out_op<<all_mtf[i].site[j].header<<"\t"
                <<all_mtf[i].site[j].strand<<"\t"<<all_mtf[i].site[j].begin
                <<"\t"<<all_mtf[i].site[j].end<<endl;
      }
      
      f_out_op<<endl<<endl;
    }
    
    f_out_op<<endl;
    
    f_out_op<<"Time "<<difftime(end_t, begin_t)<<" secs.\n";
    f_out_op.close();
    
    cout<<"Finished generating the motif file "<<tmp_out<<".\n";
    cout<<"There're altogether "<<num_mtf<<" motifs output.\n";
    cout<<"Thank you for using ProSampler!\n";
  }
  
  
  /*************************************************************************************************
   *
   * Function get_spic used to get SPIC format output file from motif
   *
   ************************************************************************************************/
  
  void get_spic(vector<mtf>& all_mtf, vector<spic>& all_spic)
  {
    spic tmp_spic;
    
    for(int i=0; i<all_mtf.size(); i++)
    {
      tmp_spic.pfm=all_mtf[i].mat;
      tmp_spic.pssm=all_mtf[i].mat;
      
      for(int j=0; j<4; j++)
      {
        for(int k=0; k<all_mtf[i].alength; k++)
        {
          tmp_spic.pfm[j][k]*=all_mtf[i].nsites;
          
          if(tmp_spic.pssm[j][k]==0)
          {
            tmp_spic.pssm[j][k]=0.00001;
          }
          
          tmp_spic.pssm[j][k]/=nt_freq[j];
          tmp_spic.pssm[j][k]=log(tmp_spic.pssm[j][k]);
        }
      }
      
      tmp_spic.ic.clear();
      
      for(int j=0; j<all_mtf[i].alength; j++)
      {
        float tmp_sum=0;
        
        for(int k=0; k<4; k++)
        {
          tmp_sum+=all_mtf[i].mat[k][j]*tmp_spic.pssm[k][j];
        }
        
        tmp_spic.ic.push_back(tmp_sum);
      }
      
      all_spic.push_back(tmp_spic);
    }
    
    cout<<"Finished transforming motifs into SPIC format.\n";
  }
  
  
  /*************************************************************************************************
   *
   * Function print_spic used to print the result into one file
   *
   ************************************************************************************************/
  
  void print_spic(char *in_file, string& bg_file, string& out_file, vector<spic>& all_spic, 
                  vector<mtf>& all_mtf)
  {
    string tmp_out=out_file;
    tmp_out+=".spic";
    ofstream f_out_op(tmp_out.c_str());
    
    if(!f_out_op)
    {
      cerr<<"Error: Can't open file "<<tmp_out.c_str()<<" for output!\n";
      exit(1);
    }
    
    f_out_op<<"# MEME 4.0.0\n"<<"# Command: ./ProSampler -i "<<in_file
            <<" -b "<<bg_file.c_str()<<" -o "<<out_file.c_str()<<" -d "<<num_deg<<" -m "<<num_mtf
            <<" -f "<<num_iter<<" -k "<<kmer_length<<" -l "<<lmer_length<<" -r "
            <<redundant_thr<<" -p "<<str_flag<<" -t "<<thr1<< " -w " << thr3 << " -c "<<hd_thr<<" -z "
            <<thr2<<" -h "<<help_flag<< " -P " << padj_cutoff << "\n";
    f_out_op<<endl<<"# Begin: "<<begin_pkg<<endl;
    f_out_op<<"#   End: "<<end_pkg<<endl<<endl;
    f_out_op<<"ProSampler Version 1.0.0\n"<<endl;
    f_out_op<<"ALPHABET= ACGT\n"<<endl;
    
    if(str_flag==1)
    {
      f_out_op<<"Strands: +\n"<<endl;
    }
    else
    {
      f_out_op<<"Strands: + -\n"<<endl;
    }
    
    f_out_op<<"Background letter frequencies (from dataset):\n";
    
    if(num_mtf>all_spic.size())
    {
      num_mtf=all_spic.size();
    }
    
    for(int i=0; i<num_mtf; i++)
    {
      f_out_op<<"MOTIF "<<all_mtf[i].deg<<" "<<all_mtf[i].rev_deg<<" ProSampler\n"<<endl;
      
      for(int j=0; j<4; j++)
      {
        f_out_op<<alphabet[j]<<"\t";
        
        for(int k=0; k<all_mtf[i].alength-1; k++)
        {
          f_out_op<<setiosflags(ios::fixed)<<setprecision(PRECISE);
          f_out_op<<setw(6)<<all_spic[i].pssm[j][k]<<"\t";
        }
        f_out_op<<setiosflags(ios::fixed)<<setprecision(PRECISE);
        f_out_op<<setw(6)<<all_spic[i].pssm[j][all_mtf[i].alength-1]<<endl;
      }
      
      f_out_op<<"I"<<"\t";
      f_out_op.setf(ios::fixed);
      
      for(int k=0; k<all_mtf[i].alength-1; k++)
      {
        f_out_op<<fixed<<setprecision(PRECISE);
        f_out_op<<setw(6)<<all_spic[i].ic[k]<<"\t";
      }
      f_out_op<<fixed<<setprecision(PRECISE);
      f_out_op<<setw(6)<<all_spic[i].ic[all_mtf[i].alength-1]<<endl;
      
      f_out_op.unsetf(ios::fixed);
      
      for(int j=0; j<4; j++)
      {
        char tmp_char;
        tmp_char=tolower(alphabet[j]);
        f_out_op<<tmp_char<<"\t";
        for(int k=0; k<all_mtf[i].alength-1; k++)
        {
          f_out_op<<int(all_spic[i].pfm[j][k])<<"\t";
        }
        f_out_op<<int(all_spic[i].pfm[j][all_mtf[i].alength-1])<<endl;
      }
      
      f_out_op<<endl;
    }
    
    f_out_op<<endl;
    
    f_out_op<<"Time "<<difftime(end_t, begin_t)<<" secs.\n";
    f_out_op.close();
    
    cout<<"Finished generating the motif file "<<tmp_out<<".\n";
    cout<<"There're altogether "<<num_mtf<<" motifs output.\n";
    cout<<"Thank you for using ProSampler!\n";
  }
  
  
  /*************************************************************************************************
   *
   * Main function
   *
   ************************************************************************************************/
  
  int main(int argc, char **argv)
  {
    begin_t = time(0);
    strftime(begin_pkg, sizeof(begin_pkg), 
             "%b %d %Y %a %X %Z", localtime(&begin_t));
    
    parse_opt(argc, argv);
    
    if(help_flag == 1)
    {
      usage();
    }
    
    r_alphabet['A'] = 0;
    r_alphabet['C'] = 1;
    r_alphabet['G'] = 2;
    r_alphabet['T'] = 3;
    
    sequence seq_nondeg;	// sequence with degenerate positions
    load_data(argv[f_in_id], seq_nondeg);
    
    // for (int i = 0; i < seq_nondeg.strand.size(); i ++)
    // {
    //   cout << seq_nondeg.strand[i] << "\n\n";
    // }
    // 
    // exit(0);
    
    // int markov_flag;
    string bg_file = "order Markov chain model";
    string OutFile = argv[f_in_id];
    // struct kmer_package kmer_pkg; // The probability model to save the information of background
    bg_file = "3-rd order Markov chain model";
    cout << "Using 3-rd order Markov Chain model as background.\n";
    // markov_flag = 3;
    struct kmer_package kmer_pkg = markov(seq_nondeg);
    
    if(f_out_id != -1)
    {
      OutFile = argv[f_out_id];
    }
    
    
    // for (int i = 0; i < bg_nondeg.strand.size(); i ++)
    // {
    //   cout << bg_nondeg.strand[i] << "\n\n";
    // }
    
    // for (map_string::iterator it = kmer_pkg.AA.begin(); it != kmer_pkg.AA.end(); it ++)
    // {
    //   cout << it -> first << "\t" << it -> second << "\n";
    // }
    
    // for (map_string::iterator it = kmer_pkg.B.begin(); it != kmer_pkg.B.end(); it ++)
    // {
    //   cout << it -> first << "\t" << it -> second << "\n";
    // }
    
    // for (map_string::iterator it = kmer_pkg.D.begin(); it != kmer_pkg.D.end(); it ++)
    // {
    //   cout << it -> first << "\t" << it -> second << "\n";
    // }
    // 
    // exit(0);
    
    de_lower(seq_nondeg);
    de_other(seq_nondeg);
    nt_stat(seq_nondeg, nt_freq);
    sequence seq_final;
    // sequence bg_final;
    seq_generate(seq_nondeg, seq_final, str_flag);
    // seq_generate(bg_nondeg, bg_final, str_flag);
    free_seq(seq_nondeg);
    // free_seq(bg_nondeg);

        
    // for (int i = 0; i < seq_final.strand.size(); i ++)
    // {
    //   cout << seq_final.strand[i] << "\n\n";
    // }
    
    // for (int i = 0; i < bg_final.strand.size(); i ++)
    // {
    //   cout << bg_final.strand[i] << "\n\n";
    // }
    // 
    // exit(0);
    
    map<string, int_vector> kmer_site;
    kmer_count(kmer_site, seq_final, lmer_length);
    cout << "Identified " << kmer_site.size() << " " << 
      kmer_length << "-mer in experiment sequences.\n";
    
    
    // for (map<string, int_vector>::iterator it = kmer_site.begin(); it != kmer_site.end(); it ++)
    // {
    //   cout << it -> first << "\t" << (it -> second)[0] << "\t" << (it -> second).size() << "\n";
    //   // break;
    // }
    
    // cout << "CAATGGCA " << kmer_site["CAATGGCA"].size() << "\n";
    
    map<string, float> kmer_bg;	// hold the k-mers in background sequences
    
    int nkmer = 0;
    
    for (map<string, int_vector>::iterator it = kmer_site.begin(); it != kmer_site.end(); it ++)
    {
      nkmer += (it -> second).size(); // Add the occurrence of this k-mer
    }
    
    kmer_in_bg(kmer_bg, kmer_pkg, kmer_site, nkmer);
    cout << "Identified " << kmer_site.size() << " " << 
      kmer_length << "-mer in background sequences.\n";
    
    
    // for (map<string, int_vector>::iterator it = kmer_bg.begin(); it != kmer_bg.end(); it ++)
    // {
    //   cout << it -> first << "\t" << (it -> second)[0] << "\t" << (it -> second).size() << "\n";
    //   // break;
    // }
    
    // cout << "CAATGGCA " << kmer_bg["CAATGGCA"].size() << "\n";
    
    // exit(0);
    
    kmer_set major_set, minor_set;	// k-mer sets
    map<string, int> fa_flag;	
    // The flag to represent the status of k-mers in ChIP-seq data
    map<string, int> bg_flag;	// Status in background data
    
    for(map<string, int_vector>::iterator it = kmer_site.begin(); it != kmer_site.end(); it ++)
    {
      fa_flag[it -> first] = 1;
    }
    
    // Combine k-mers in positive and negative strands
    if(str_flag != 1)
    {
      comb_kmer(kmer_site, fa_flag, seq_final.strand, 
                kmer_bg, bg_flag);
    }
    
    // cout << "CAATGGCA " << kmer_site["CAATGGCA"].size() << "\n";
    // cout << "CAATGGCA " << kmer_bg["CAATGGCA"].size() << "\n";
    // exit(0);
    int data_size = 0;
    for (int i = 0; i < seq_final.strand.size(); i ++) {
      data_size += seq_final.strand[i].length();
    }
    cout << "The total size of the data is " << data_size << ".\n";
    if (data_size <= correction_cutoff) {
      cout << "We will not conduct correction for p-values ...\n";
      correction_method = 0;
    }
    
    choose_kmer(kmer_site, kmer_bg, major_set, correction_method, 
                seq_final.name.size(), fa_flag, minor_set, thr1, thr3);
    // select the significant and sub-significant k-mers
    
    kmer_sort(major_set); // sort all significant k-mers
    
    if (correction_method == 2) {
      BH_correct(major_set, minor_set, thr1, thr3);
    }
    
    // cout << major_set.kmer[10] << " " << major_set.occurr[10] << " " << major_set.cover1[10] << " " <<
    //   major_set.cover2[10] << " " << major_set.score[10] << " " << major_set.site[10].size() << "\n";
    // 
    // exit(0);
    
    if (correction_method != 2) {
      kmer_sort(minor_set);
    }
    
    kmer_site.clear();
    // kmer_bg.clear();
    map<string, int_vector>().swap(kmer_site);
    map<string, float>().swap(kmer_bg);
    // free_seq(bg_nondeg);
    // free_seq(bg_final);
    
    // for (int i = 0; i < major_set.kmer.size(); i ++)
    // {
    //   cout << major_set.kmer[i] << " " << major_set.occurr[i] << " " << major_set.cover1[i] << " " <<
    //     major_set.cover2[i] << " " << major_set.score[i] << " " << major_set.site[i].size() << "\n";
    // }
    // 
    // exit(0);
    
    // cout << major_set.kmer[10] << " " << major_set.occurr[10] << " " << major_set.cover1[10] << " " <<
    //   major_set.cover2[10] << " " << major_set.score[10] << " " << major_set.site[10].size() << "\n";
    // 
    // exit(0);
    
    // for (int i = 0; i < minor_set.kmer.size(); i ++)
    // {
    //   cout << minor_set.kmer[i] << " " << minor_set.occurr[i] << " " << minor_set.cover1[i] << " " <<
    //     minor_set.cover2[i] << " " << minor_set.score[i] << " " << minor_set.site[i].size() << "\n";
    // }
    // 
    // exit(0);
    
    top_num = combine_major(major_set, minor_set);
    free_kmer(minor_set);
    vector<kmer_all> kmer_final;
    fill_kmer(major_set, kmer_final, seq_final.name.size(), seq_final);
    
    // for (int i = 0; i < kmer_final.size(); i ++)
    // {
    //   cout << kmer_final[i].kmer << " " << kmer_final[i].occurr << " " << kmer_final[i].cover1 <<
    //     " " << kmer_final[i].cover2 << " " << kmer_final[i].score << "\n";
    // }
    // 
    // exit(0);
    
    free_kmer(major_set);
    vector<psm> psm_set;
    kmer2psm(kmer_final, psm_set, top_num);
    
    // 	for (int i = 0; i < psm_set.size(); i ++)
    // 	{
    // 	  for (uint j = 0; j < psm_set[i].set.size(); j ++)
    // 	  {
    // 	    cout << psm_set[i].set[j] << " ";
    // 	  }
    // 	  
    // 	  cout << "\n";
    // 
    // 	  // for (int j = 0; j < psm_set[i].mat.size(); j ++)
    // 	  // {
    // 	  //   for (int k = 0; k < psm_set[i].mat[j].size(); k ++)
    // 	  //   {
    // 	  //     cout << psm_set[i].mat[j][k] << " ";
    // 	  //   }
    // 	  // 
    // 	  //   cout << "\n";
    // 	  // }
    // // 
    // // 	  cout << "\n\n\n\n";
    // 	}
    // 	
    // 	exit(0);
    
    fill_psm(psm_set, kmer_final, kmer_length); // compute the matrices for all PSMs
    
    // for (int i = 0; i < psm_set.size(); i ++)
    // {
    //   for (int j = 0; j < psm_set[i].mat.size(); j ++)
    //   {
    //     for (int k = 0; k < psm_set[i].mat[j].size(); k ++)
    //     {
    //       cout << psm_set[i].mat[j][k] << " ";
    //     }
    // 
    //     cout << "\n";
    //   }
    // 
    //   cout << "\n\n\n\n";
    // }
    // 
    // exit(0);
    
    int seq_num=seq_final.name.size();
    cover_psm(psm_set, kmer_final);
    test_psm(psm_set, sum_exp, sum_bg);
    comp_cons(psm_set);
    psm_sort(psm_set);
    vector<int_vector> pre_out;
    int num_pre;
    
    // for (int i = 0; i < psm_set.size(); i ++)
    // {
    //   cout << psm_set[i].occurr << " " << psm_set[i].cons << " " << psm_set[i].cover1 << " " << 
    //     psm_set[i].cover2 << "\n";
    // }
    // 
    // exit(0);
    
    if(num_mtf != -1)
    {
      num_pre=MULTIPLE*num_mtf;
    }
    else
    {
      num_pre=psm_set.size();
    }
    
    vector<int_vector> psm_nbr;
    build_graph(psm_set, psm_nbr);
    gibbs_sam(pre_out, psm_set, kmer_final,
              num_pre, num_iter, psm_nbr, seq_num);
    
    // for (int i = 0; i < pre_out.size(); i ++)
    // {
    //   for (int j = 0; j < pre_out[i].size(); j ++)
    //   {
    //     cout << pre_out[i][j] << " ";
    //   }
    //   
    //   cout << "\n";
    // }
    // 
    // exit(0);
    
    psm_set.clear();
    vector<psm>().swap(psm_set);
    psm_nbr.clear();
    vector<int_vector>().swap(psm_nbr);
    vector<pre_mtf> all_pre;
    get_pre(all_pre, pre_out, kmer_final, seq_num);
    
    // for (int i = 0; i < all_pre.size(); i ++)
    // {
    //   for (int j = 0; j < all_pre[i].set.size(); j ++)
    //   {
    //     cout << j << " ";
    //   }
    //   cout << "\n";
    // }
    // 
    // exit(0);
    
    pre_out.clear();
    vector<int_vector>().swap(pre_out);
    refine_mtf(all_pre);
    
    // exit(0);
    
    for (map<string, float>::iterator it = kmer_bg.begin(); it != kmer_bg.end(); it ++) {
      (it -> second) /= nkmer;
    }
    compute_pvalue(all_pre, kmer_final, kmer_bg); // sift preliminary motifs using chi square p-values
    vector<kmer_all>().swap(kmer_final);
    
    // exit(0);
    
    kmer_final.clear();
    make_wild(wild_card, reverse_wild);
    vector<mtf> all_mtf;
    pre2mtf(all_pre, all_mtf, seq_final);
    all_pre.clear();
    vector<pre_mtf>().swap(all_pre);
    sort_mtf(all_mtf);
    sift_mtf(all_mtf);
    
    end_t = time(0);
    strftime(end_pkg, sizeof(end_pkg), 
             "%b %d %Y %a %X %Z", localtime(&end_t));
    print_mtf(argv[f_in_id], bg_file, OutFile, all_mtf);
    print_site(argv[f_in_id], bg_file, OutFile, all_mtf);
    vector<spic> all_spic;
    all_spic.reserve(all_mtf.size());
    get_spic(all_mtf, all_spic);
    print_spic(argv[f_in_id], bg_file, OutFile, all_spic, all_mtf);
    
    return 0;
  }
  