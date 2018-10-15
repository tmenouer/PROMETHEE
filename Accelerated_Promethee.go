/*
 *
 * Accelerated PROMETHEE Algorithm based on K-means
 *
 * For the DTW distance, we use https://github.com/mjanda/go-dtw
 * This is the only one dependency
 *
 * Christophe Cerin and Tarek Menouer - Oct 2018
 */

package main

import (
	"container/list"
	//"encoding/json"
	"fmt"
	//"github.com/Arafatk/glot"
	"flag"
	"math"
	"math/rand"
	//"reflect"
	"sync"
	"time"
	"sort"
	"go-dtw"
)

type info struct {
	cpu      int
	memory   int
	f_cpu    int
	f_memory int
	id_value int
}

type node struct {
	center_cpu    float64
	center_memory float64
	size          int
	membres       []info
}

type Struct_PROMETHEE struct {
	D []float64
	P float64
}

type ByPromethee_value []map[int]float64

func (a ByPromethee_value) Len() int           { return len(a) }
func (a ByPromethee_value) Swap(i, j int)      { a[i], a[j] = a[j], a[i] }
func (a ByPromethee_value) Less(i, j int) bool { return a[i][2*SIZE+1] > a[j][2*SIZE+1] }
var SIZE int
func Kmeans(MyTab []map[int]int, k int, taille int, nb_elements int) []map[int]int {

	var init, cmpt, id, somme_cpu, somme_memory, somme_elements, nb_elements_cluster int = 0, 0, 0, 0, 0, 0, 0
	var stop_init, stop_change bool = false, true
	var distance_old, distance float64 = 0, 0

	Resultat := make([]node, k)
	Old_Resultat := make([]node, k)
	TableInit := make([]int, k)

	for i := 0; i < len(TableInit); i++ {
		TableInit[i] = -1
	}

	//Initialiser les centres
	for {
		init = rand.Intn(taille)
		stop_init = false
		for j := 0; j < len(TableInit); j++ {
			if TableInit[j] == init {
				stop_init = true
			}
		}

		if stop_init == false {

			Resultat[cmpt].center_cpu = float64(MyTab[init][0])
			Resultat[cmpt].center_memory = float64(MyTab[init][2])
			Resultat[cmpt].size = 0
			Resultat[cmpt].membres = make([]info, taille)

			Old_Resultat[cmpt].center_cpu = float64(MyTab[init][0])
			Old_Resultat[cmpt].center_memory = float64(MyTab[init][2])
			Old_Resultat[cmpt].size = 0
			Old_Resultat[cmpt].membres = make([]info, taille)

			cmpt++
		}

		if cmpt == k {
			break
		}
	}

	//Boucle principale
	for {
		//Calculer la distance
		for i := 0; i < taille; i++ {
			distance_old = math.Abs(float64(MyTab[i][0])-Resultat[0].center_cpu) + math.Abs(float64(MyTab[i][2])-Resultat[0].center_memory)
			id = 0
			for j := 1; j < k; j++ {
				distance = math.Abs(float64(MyTab[i][0])-Resultat[j].center_cpu) + math.Abs(float64(MyTab[i][2])-Resultat[j].center_memory)
				if distance < distance_old {
					distance_old = distance
					id = j
				}
			}

			Resultat[id].membres[Resultat[id].size].cpu = MyTab[i][0]
			Resultat[id].membres[Resultat[id].size].f_cpu = MyTab[i][1]
			Resultat[id].membres[Resultat[id].size].memory = MyTab[i][2]
			Resultat[id].membres[Resultat[id].size].f_memory = MyTab[i][3]
			Resultat[id].membres[Resultat[id].size].id_value = MyTab[i][4]
			Resultat[id].size = Resultat[id].size + 1
		}

		//Calcul les nouveaux centres
		for i := 0; i < k; i++ {
			somme_cpu = 0
			somme_memory = 0
			for j := 0; j < Resultat[i].size; j++ {
				somme_cpu += Resultat[i].membres[j].cpu
				somme_memory += Resultat[i].membres[j].memory
			}
			if Resultat[i].size != 0 {
				Resultat[i].center_cpu = float64(somme_cpu / Resultat[i].size)
				Resultat[i].center_memory = float64(somme_memory / Resultat[i].size)
			} else {
				Resultat[i].center_cpu = 0
				Resultat[i].center_memory = 0
			}
		}

		//Voir s'il y a un changement
		stop_change = true
		for i := 0; i < k; i++ {
			if (Resultat[i].center_cpu != Old_Resultat[i].center_cpu) || (Resultat[i].center_memory != Old_Resultat[i].center_memory) {
				stop_change = false
			} else {
				if Resultat[i].size != Old_Resultat[i].size {
					stop_change = false
				} else {
					for j := 0; j < Resultat[i].size; j++ {
						if (Resultat[i].membres[j].cpu != Old_Resultat[i].membres[j].cpu) || (Resultat[i].membres[j].memory != Old_Resultat[i].membres[j].memory) {
							stop_change = false
						}
					}
				}
			}
		}

		//Oui changeme nt
		if stop_change == false {
			for i := 0; i < k; i++ {

				Old_Resultat[i].center_cpu = Resultat[i].center_cpu
				Old_Resultat[i].center_memory = Resultat[i].center_memory
				Old_Resultat[i].size = Resultat[i].size
				Old_Resultat[i].membres = nil
				Old_Resultat[i].membres = make([]info, taille)

				for j := 0; j < Resultat[i].size; j++ {
					Old_Resultat[i].membres[j].cpu = Resultat[i].membres[j].cpu
					Old_Resultat[i].membres[j].memory = Resultat[i].membres[j].memory
				}
				Resultat[i].size = 0
				Resultat[i].membres = nil
				Resultat[i].membres = make([]info, taille)
			}
		} else { //Non fin
			break
		}
	}

	//Calcul le nombre total d'elements
	for i := 0; i < k; i++ {
		if Resultat[i].size > nb_elements {
			nb_elements_cluster = nb_elements
		} else {
			nb_elements_cluster = Resultat[i].size
		}
		somme_elements += nb_elements_cluster
	}

	//insertion de nb_elements de chaque cluster
	Return_MyTab := make([]map[int]int, somme_elements)
	for i := range Return_MyTab {
		Return_MyTab[i] = make(map[int]int)
	}
	cmpt = 0
	for i := 0; i < k; i++ {
		if Resultat[i].size > nb_elements {
			nb_elements_cluster = nb_elements
		} else {
			nb_elements_cluster = Resultat[i].size
		}
		for j := 0; j < nb_elements_cluster; j++ {
			Return_MyTab[cmpt][0] = Resultat[i].membres[j].cpu
			Return_MyTab[cmpt][1] = Resultat[i].membres[j].f_cpu
			Return_MyTab[cmpt][2] = Resultat[i].membres[j].memory
			Return_MyTab[cmpt][3] = Resultat[i].membres[j].f_memory
			Return_MyTab[cmpt][4] = Resultat[i].membres[j].id_value
			cmpt++
		}
	}
	
	return Return_MyTab
}


type Value float64

//
//  Helper sorting functions.
//

// is v < w ?
func less(v int, w int) bool {
	if v == w {
		return false
	} else // optimization when reference equal
	{
		return (v < w)
	}
}



// does v == w ?
func eq(v int, w int) bool {
	if v == w {
		return true
	} // optimization when reference equal
	return false
}

// exchange a[i] and a[j]
func exch(a []int, i int, j int) {
	swap := a[i]
	a[i] = a[j]
	a[j] = swap
}

//
//
// 3way Quicksort from Robert Sedgewick
// https://www.cs.princeton.edu/~rs/talks/QuicksortIsOptimal.pdf
//

func _3way(a []int, l int, r int) {
	if r <= l {
		return
	}
	v := a[r]
	i := l - 1
	j := r
	p := l - 1
	q := r
	for {
		i++
		for less(a[i], v) {
			i++
		}
		j--
		for less(v, a[j]) {
			if j == l {
				break
			} else {
				j--
			}
		}
		if i >= j {
			break
		}
		exch(a, i, j)
		if eq(a[i], v) {
			p++
			exch(a, p, i)
		}
		if eq(v, a[j]) {
			q--
			exch(a, q, j)
		}
	}
	exch(a, i, r)
	j = i - 1
	i = i + 1
	for k := l; k < p; k, j = k+1, j-1 {
		exch(a, k, j)
	}
	for k := r - 1; k > q; k, i = k-1, i+1 {
		exch(a, k, i)
	}
	_3way(a, l, j)
	_3way(a, i, r)
}

//
// 3way quicksort for arrays of dictionaries/maps
// The dimension for the sort is given as a parameter (dim)
//

// is v < w ?
func lessD(v map[int]int, w map[int]int, dim int) bool {
	if v[dim] == w[dim] {
		return false
	} else // optimization when reference equal
	{
		return (v[dim] < w[dim])
	}
}

// is v > w ?
func gtD(v map[int]int, w map[int]int, dim int) bool {
	if v[dim] == w[dim] {
		return false
	} else // optimization when reference equal
	{
		return (v[dim] > w[dim])
	}
}

// does v == w ?
func eqD(v map[int]int, w map[int]int, dim int) bool {
	if v[dim] == w[dim] {
		return true
	} // optimization when reference equal
	return false
}

// exchange a[i] and a[j]
func exchD(a []map[int]int, i int, j int, dim int) {
	swap := map[int]int{}
	for key, _ := range a[i] {
		swap[key] = a[i][key]
	}
	for key, _ := range a[i] {
		a[i][key] = a[j][key]
	}
	for key, _ := range a[i] {
		a[j][key] = swap[key]
	}
}

func InsertionSortD(v []map[int]int, dim int) {
	for j := 1; j < len(v); j++ {
		// Invariant: v[:j] contains the same elements as
		// the original slice v[:j], but in sorted order.
		key := v[j]
		i := j - 1
		for i >= 0 && gtD(v[i], key, dim) {
			v[i+1] = v[i]
			i--
		}
		v[i+1] = key
	}
}

func _3wayD(a []map[int]int, l int, r int, dim int) {
	//if r <= l {
	//	return
	//}
	if r-l < 4096 {
		InsertionSortD(a[l:r+1], dim)
		return
	}
	v := make(map[int]int)
	for key, value := range a[r] {
		// fmt.Println("Key:", key, "Value:", value)
		v[key] = value
	}
	i := l - 1
	j := r
	p := l - 1
	q := r
	for {
		i++
		for lessD(a[i], v, dim) {
			i++
		}
		j--
		for lessD(v, a[j], dim) {
			if j == l {
				break
			} else {
				j--
			}
		}
		if i >= j {
			break
		}
		exchD(a, i, j, dim)
		if eqD(a[i], v, dim) {
			p++
			exchD(a, p, i, dim)
		}
		if eqD(v, a[j], dim) {
			q--
			exchD(a, q, j, dim)
		}
	}
	exchD(a, i, r, dim)
	j = i - 1
	i = i + 1
	for k := l; k < p; k, j = k+1, j-1 {
		exchD(a, k, j, dim)
	}
	for k := r - 1; k > q; k, i = k-1, i+1 {
		exchD(a, k, i, dim)
	}
	_3wayD(a, l, j, dim)
	_3wayD(a, i, r, dim)
}

func _3wayD_PAR(a []map[int]int, l int, r int, dim int) {
	if r <= l {
		return
	}
	v := make(map[int]int)
	for key, value := range a[r] {
		// fmt.Println("Key:", key, "Value:", value)
		v[key] = value
	}
	i := l - 1
	j := r
	p := l - 1
	q := r
	for {
		i++
		for lessD(a[i], v, dim) {
			i++
		}
		j--
		for lessD(v, a[j], dim) {
			if j == l {
				break
			} else {
				j--
			}
		}
		if i >= j {
			break
		}
		exchD(a, i, j, dim)
		if eqD(a[i], v, dim) {
			p++
			exchD(a, p, i, dim)
		}
		if eqD(v, a[j], dim) {
			q--
			exchD(a, q, j, dim)
		}
	}
	exchD(a, i, r, dim)
	j = i - 1
	i = i + 1
	for k := l; k < p; k, j = k+1, j-1 {
		exchD(a, k, j, dim)
	}
	for k := r - 1; k > q; k, i = k-1, i+1 {
		exchD(a, k, i, dim)
	}

	var wg sync.WaitGroup
	wg.Add(2)
	go func() { _3wayD_PAR(a, l, j, dim); wg.Done() }()
	go func() { _3wayD_PAR(a, i, r, dim); wg.Done() }()
	wg.Wait()
}


func PROMETHEE(MyTab []map[int]int, nodeCount int, nbDimension int) *list.List {
	
	Promethee_MyTab := make([]map[int]float64, nodeCount)
			
	for i := range Promethee_MyTab {
		Promethee_MyTab[i] = make(map[int]float64)
		
		for j:=0;j<((2*nbDimension)+1);j++ {
			Promethee_MyTab[i][j]=float64(MyTab[i][j])
		}
		Promethee_MyTab[i][(2*nbDimension)+1] = 0
	}
	
	
	Matrix_Nodes := make([][]Struct_PROMETHEE, nodeCount)

	for i := range Matrix_Nodes {
		Matrix_Nodes[i] = make([]Struct_PROMETHEE, nodeCount)
	}

	var fp, fn float64

	for i := 0; i < nodeCount; i++ {
		for j := 0; j < nodeCount; j++ {
			Matrix_Nodes[i][j].D = make([]float64, nbDimension)
			Matrix_Nodes[i][j].P = 0
			
			for k := 0; k < nbDimension; k++ {
				Matrix_Nodes[i][j].D[k] = Promethee_MyTab[i][k*2] - Promethee_MyTab[j][k*2]
				
				if Matrix_Nodes[i][j].D[k] > 0 {
					Matrix_Nodes[i][j].P += 1
				}
			}
			//~ Matrix_Nodes[i][j].D[0] = Promethee_MyTab[i][0] - Promethee_MyTab[j][0]
			//~ Matrix_Nodes[i][j].D[1] = Promethee_MyTab[i][2] - Promethee_MyTab[j][2]
		}
	}

	for i := 0; i < nodeCount; i++ {
		fp = 0
		fn = 0
		for j := 0; j < nodeCount; j++ {
			fp = fp + Matrix_Nodes[i][j].P
		}
		fp = fp / float64(nodeCount-1)

		for j := 0; j < nodeCount; j++ {
			fn = fn + Matrix_Nodes[j][i].P
		}

		fn = fn / float64(nodeCount-1)

		//l[i].Promethee_value = fp - fn
		Promethee_MyTab[i][2*nbDimension+1] = fp - fn
	}
	
	sort.Sort(ByPromethee_value(Promethee_MyTab))
	
	R_Promethee_MyTab := make([]map[int]int, nodeCount)
			
	for i := range R_Promethee_MyTab {
		R_Promethee_MyTab[i] = make(map[int]int)
		
		for j:=0;j<(2*nbDimension)+1;j++ {
			R_Promethee_MyTab[i][j]=int(Promethee_MyTab[i][j])
		}
		
	}

	l := list.New()

	for i := 0; i < nodeCount; i++ {
		l.PushBack(R_Promethee_MyTab[i])
	}

	return l
	
}

func Maximization_objective(MyTab []map[int]int) []map[int]int {

	if MyTab[0][1] == 1 {

		//If function 1 is maximization
		x := make([]int, 5)
		for i := 0; i < len(MyTab)/2; i++ {
			for j := 0; j < 5; j++ {
				x[j] = MyTab[i][j]
			}

			for j := 0; j < 5; j++ {
				MyTab[i][j] = MyTab[len(MyTab)-1-i][j]
			}

			for j := 0; j < 5; j++ {
				MyTab[len(MyTab)-1-i][j] = x[j]
			}
		}
	}
	return MyTab
}

func present(x int, y int, l *list.List) bool {
	e := l.Front()
	for index1 := 0; index1 < l.Len(); index1++ {
		if x == e.Value.(map[int]int)[0] && y == e.Value.(map[int]int)[2] {
			return true
		}
		e = e.Next()
	}
	return false
}

func JaccardIndex(l_exact []*list.List, l_proj []*list.List, taille int, n int, m int) {
	my_mod := int(taille / n)
	fmt.Println("Size of buffers:", taille, "Modulo:", my_mod, "Min:", m)
	result := make([]map[string]float64, my_mod)
	for i := range result {
		result[i] = map[string]float64{"jaccard": 0.0}
	}
	l_index := list.New()
	l_index1 := list.New()
	for index := 0; index < taille; index++ {
		l_index = l_proj[index]
		l_index1 = l_exact[index]
		e := l_index.Front()
		union := l_index.Len()
		intersect := 0
		for index1 := 0; index1 < l_index.Len(); index1++ {
			//fmt.Println(index,index1)
			//for key, value := range e.Value.(map[int]int) {
			//    fmt.Println("Key:", key, "Value:", value)
			//}
			//fmt.Println("Val0: ",e.Value.(map[int]int)[0])
			//fmt.Println("Val1: ",e.Value.(map[int]int)[1])
			//fmt.Println("Val2: ",e.Value.(map[int]int)[2])
			//fmt.Println("Val3: ",e.Value.(map[int]int)[3])
			//fmt.Println("Val4: ",e.Value.(map[int]int)[4])
			if present(e.Value.(map[int]int)[0], e.Value.(map[int]int)[2], l_index1) {
				intersect += 1
			} else {
				union += 1
			}
			e = e.Next()
		}
		result[index%my_mod]["jaccard"] = (result[index%my_mod]["jaccard"] + float64(intersect)/float64(union)) / 2.0
		//fmt.Println(index,"Union",union,"Interset",intersect,"Jaccard Index",result[index%my_mod]["jaccard"])
	}
	// Print results
	for i := range result {
		fmt.Println("Jaccard Index for size", int(float64(m)*math.Pow(float64(2), float64(i))), ":", result[i]["jaccard"])
	}
}

func distance(x int, y int) float64 {
	// Euclidien distance in 2D space
	s := 0.0
	d := float64(x)
	s += d * d
	d = float64(y)
	s += d * d
	return math.Sqrt(s)
}

func JaccardSimilarityCoeff(l_exact []*list.List, l_proj []*list.List, taille int, n int, m int) {
	// We compute the Jaccard Similarity coefficients. We assume that the 2 lists/buffers
	// have equal sizes.
	my_mod := int(taille / n)
	fmt.Println("Size of buffers:", taille, "Modulo:", my_mod, "Min:", m)
	result := make([]map[string]float64, my_mod)
	for i := range result {
		result[i] = map[string]float64{"min": 0.0, "max": 0.0}
	}
	l_index1 := list.New()
	l_index2 := list.New()
	for index := 0; index < taille; index++ {
		l_index1 = l_proj[index]
		l_index2 = l_exact[index]
		e1 := l_index1.Front()
		e2 := l_index2.Front()
		d1 := distance(e1.Value.(map[int]int)[0], e1.Value.(map[int]int)[2])
		d2 := distance(e2.Value.(map[int]int)[0], e2.Value.(map[int]int)[2])
		if d1 > d2 {
			result[index%my_mod]["min"] = result[index%my_mod]["min"] + d2
			result[index%my_mod]["max"] = result[index%my_mod]["max"] + d1
		} else {
			result[index%my_mod]["min"] = result[index%my_mod]["min"] + d1
			result[index%my_mod]["max"] = result[index%my_mod]["max"] + d2
		}
		e1 = e1.Next()
		e2 = e2.Next()
	}
	// Print results
	for i := range result {
		fmt.Println("Jaccard Similarity Coefficient for size", int(float64(m)*math.Pow(float64(2), float64(i))), ":", float64(result[i]["min"]/result[i]["max"]))
	}
}

func main() {
	// Controlling the command line
	flag.Parse()
	if flag.NArg() != 1 {
		fmt.Println("Bag number of arguments\n")
		fmt.Println("\tFlag: kmeans to use the kmeans method")
		fmt.Println("\tFlag: default to use the default Exact Promethee method")
		fmt.Println("\tFlag: analysis to compute similarity metrics")
		return
	}
	fmt.Println("Arguments: ", flag.Args())
	if flag.Args()[0] == "kmeans" {
		fmt.Println("Kmeans method selected")
	} else {
		if flag.Args()[0] == "default" {
			fmt.Println("Exact Promethee method selected")
		} else {
			if flag.Args()[0] == "analysis" {
				fmt.Println("We compute the similarity metrics of performance")
			} else {
				fmt.Println("Bag argument\n\tFlag: kmeans to use the kmeans method")
				fmt.Println("\tFlag: default to use the default Exact Promethee method")
				fmt.Println("\tFlag: analysis to compute similary metrics")
				return
			}
		}
	}

	//
	// Parameters of the experiment
	//
	//~ powof2 := int(math.Pow(2.0, 14.0)) // Max size ;
	powof2 := 8192 // Max size ; 16384
	MinSize := 4096                   // Min size = 4096 by default
	Dim1 := 100                   // random values between 0 and Dim1
	//~ Dim2 := 10000000                   // random values between 0 and Dim2
	const N int = 3          // Number of runs
	const NN int = int(N) + 1          // Number of runs + 1
	const Scale = 2048 //Scale between experiences
	X := ((powof2-MinSize)/Scale)+1 // Number of experiences
	SIZE = 4 //Size of multi cretiria 

	fmt.Println("Begin...")
	rand.Seed(time.Now().UTC().UnixNano())
	var result_3way [NN]map[int]float64
	var result_3way_exact [NN]map[int]float64
	var result_kung [NN]map[int]float64
	var result_kung_exact [NN]map[int]float64
	var MyTabResult []map[int]int
	// For storing the similarity metrics
	l_exact := make([]*list.List, N*X)
	l_proj := make([]*list.List, N*X)
	My_ind := 0

	for j := 0; j < int(N)+1; j++ {
		result_3way[j] = make(map[int]float64)
		result_3way_exact[j] = make(map[int]float64)
		result_kung[j] = make(map[int]float64)
		result_kung_exact[j] = make(map[int]float64)
	}

	for j := 0; j < int(N); j++ {
		for taille := MinSize; taille <= powof2; taille = taille + Scale {
			MyTab := make([]map[int]int, taille)
			l := list.New()
			for i := range MyTab {
				MyTab[i] = make(map[int]int)
				for k :=0; k<2*SIZE; k+=2 {
					MyTab[i][k] = rand.Intn(Dim1) //Value in x for function 1
					MyTab[i][k+1] = 0               // O if minimization - 1 if maximization
				}
				MyTab[i][2*SIZE] = i               //Id value
			}

			fmt.Println("Problem size ", taille)

			dimension := 0
			start := time.Now()
			if flag.Args()[0] == "default" || flag.Args()[0] == "analysis" {
				_3wayD_PAR(MyTab, 0, len(MyTab)-1, dimension)

			}
			elapsed_exact := time.Since(start)
			
			if flag.Args()[0] == "default" || flag.Args()[0] == "analysis" {	
				fmt.Printf("3way exact Quicksort parallel takes %s\n", elapsed_exact)
			}
			if flag.Args()[0] == "kmeans" || flag.Args()[0] == "analysis" {
			   /*
			    * Kmeans algorithm: Input :
			    * MyTab -- MyTab de base 
			    * K -- nombre de clusters
			    * taille -- taille de MyTab
			    * nb_elements -- nombre d'éléments à récupérer dans chaque cluster
			    *
			    * Kmeans takes taille elements and returns k*nb_elements
			    */
			    MyTabResult = Kmeans(MyTab,5,taille,50 * int(math.Log2(float64(taille))))
			    // We sort the result of the Kmeans
			    _3wayD_PAR(MyTabResult, 0, len(MyTabResult)-1, dimension)
			    
			}
			elapsed := time.Since(start)

			if flag.Args()[0] == "default" {
				fmt.Printf("3way exact parallel takes %s\n", elapsed)
			} else {
				if flag.Args()[0] == "kmeans" {
				fmt.Printf("Keans takes %s (Number of elements: %d)\n", elapsed, len(MyTabResult))
				}else
				{
					fmt.Printf("Analysis takes %s (Number of elements: %d)\n", elapsed-elapsed_exact, len(MyTabResult))	
				}
			}
			result_3way[j][taille] = float64(elapsed / time.Millisecond)
			result_3way[N][taille] = 0.0
			
			if flag.Args()[0] == "analysis" {
				result_3way_exact[j][taille] = float64(elapsed_exact / time.Millisecond)
				result_3way_exact[N][taille] = 0.0
				
				result_3way[j][taille] = float64((elapsed-elapsed_exact) / time.Millisecond)
				result_3way[N][taille] = 0.0

			}
			
			
			if flag.Args()[0] == "default" {
				MyTab = Maximization_objective(MyTab)
			} else {
			  if flag.Args()[0] != "kmeans" {
				MyTabResult = Maximization_objective(MyTabResult)
			  }
			}

			//l = KUNG(MyTab, len(MyTab), 5)
			start = time.Now()
			if flag.Args()[0] == "analysis" {
				fmt.Printf("JE SUIS AVEC %d \n",My_ind)
				l_exact[My_ind] = PROMETHEE(MyTab, len(MyTab), SIZE)
				fmt.Printf("FIN EXACT\n")
				elapsed_exact = time.Since(start)
				l_proj[My_ind] = PROMETHEE(MyTabResult, len(MyTabResult), SIZE)
				fmt.Printf("FINPROJECT\n")
				My_ind++
			} else {
				if flag.Args()[0] == "default" {
					l = PROMETHEE(MyTab, len(MyTab), SIZE)
				} else {
					l = PROMETHEE(MyTabResult, len(MyTabResult), SIZE)
				}
			}
			elapsed = time.Since(start)
			if flag.Args()[0] != "analysis" {
				fmt.Printf("Kung parallel takes %s (Number of elements in the front: %d)\n", elapsed, l.Len())
			}else
			{
				fmt.Printf("Kung parallel exact takes %s \n", elapsed_exact)
				fmt.Printf("Kung parallel project takes %s \n", elapsed - elapsed_exact)
			}
			
			if flag.Args()[0] != "analysis" {
				result_kung[j][taille] = float64(elapsed / time.Millisecond)
				result_kung[N][taille] = 0.0
			}else
			{
				result_kung_exact[j][taille] = float64(elapsed_exact / time.Millisecond)
				result_kung_exact[N][taille] = 0.0
				
				result_kung[j][taille] = float64((elapsed-elapsed_exact) / time.Millisecond)
				result_kung[N][taille] = 0.0
				
			}
		}
	}

	// We print a summary of results
	fmt.Println("-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-")
	if flag.Args()[0] == "projection" {
		fmt.Println("Mean execution time (N experiments) for Random projection")
	} else {
		if flag.Args()[0] == "default" {
			fmt.Println("Mean execution time (N experiments) for 3way")
		}
	}
	
	
	if flag.Args()[0] != "analysis" {
		for j := 0; j < int(N); j++ {
			for key, value := range result_3way[j] {
				//fmt.Println("Key:", key, "Value:", value)
				result_3way[N][int(key)] += float64(value)
			}
		}
		for key, value := range result_3way[N] {
			fmt.Println("Key:", key, "Value:", value/float64(N))
		}

		fmt.Println("Mean execution time (N experiments) for Kung")
		for j := 0; j < int(N); j++ {
			//s := 0.0
			for key, value := range result_kung[j] {
				//fmt.Println("Key:", key, "Value:", value)
				result_kung[N][int(key)] += float64(value)
			}
		}
		for key, value := range result_kung[N] {
			fmt.Println("Key:", key, "Value:", value/float64(N))
		}
	}else
	{
		for j := 0; j < int(N); j++ {
			for key, value := range result_3way[j] {
				//fmt.Println("Key:", key, "Value:", value)
				result_3way[N][int(key)] += float64(value)
			}
		}
		
		//~ for key, value := range result_3way[N] {
			//~ fmt.Println("Key:", key, "Value:", value/float64(N))
		//~ }

		//~ fmt.Println("Mean execution time (N experiments) for Kung")
		for j := 0; j < int(N); j++ {
			//s := 0.0
			for key, value := range result_kung[j] {
				//fmt.Println("Key:", key, "Value:", value)
				result_kung[N][int(key)] += float64(value)
			}
		}
		//~ for key, value := range result_kung[N] {
			//~ fmt.Println("Key:", key, "Value:", value/float64(N))
		//~ }
		
		
		for key_3way, value_3way := range result_3way[N] {
			for key, value := range result_kung[N] {
				if key_3way == key {
					fmt.Println("Proj Key:", key, "Value:", value_3way/float64(N), "value 2 ",value/float64(N), "Average Time ", value_3way/float64(N)+value/float64(N) )
				}
			}
			
		}
				
		////
		
		for j := 0; j < int(N); j++ {
			for key, value := range result_3way_exact[j] {
				//fmt.Println("Key:", key, "Value:", value)
				result_3way_exact[N][int(key)] += float64(value)
			}
		}
		//~ for key, value := range result_3way_exact[N] {
			//~ fmt.Println("Key exact:", key, "Value:", value/float64(N))
		//~ }

		//~ fmt.Println("Mean execution time (N experiments) for Kung")
		for j := 0; j < int(N); j++ {
			//s := 0.0
			for key, value := range result_kung_exact[j] {
				//fmt.Println("Key:", key, "Value:", value)
				result_kung_exact[N][int(key)] += float64(value)
			}
		}
		//~ for key, value := range result_kung_exact[N] {
			//~ fmt.Println("Key exact:", key, "Value:", value/float64(N))
		//~ }
		
		for key_3way, value_3way := range result_3way_exact[N] {
			for key, value := range result_kung_exact[N] {
				if key_3way == key {
					fmt.Println("Exact Key:", key, "Value:", value_3way/float64(N), "value 2 ",value/float64(N), "Average Time ", value_3way/float64(N)+value/float64(N) )
				}
			}
			
		}
		
	}
	
	
	if flag.Args()[0] == "analysis" {
		fmt.Println("We compute now the Jaccard indexes")
		JaccardIndex(l_exact, l_proj, N*X, N, MinSize)
		fmt.Println("We compute now the Jaccard Similarity Coefficients")
		JaccardSimilarityCoeff(l_exact, l_proj, N*X, N, MinSize)

		// prepare arrays
		a := make([]float64,N*X)
		b := make([]float64,N*X)
		//a := []float64{1, 1, 1, 2, 2, 2, 3, 3, 3, 2, 2, 4, 4, 4, 4}
		//b := []float64{1, 1, 2, 2, 3, 3, 2, 4, 4, 4}
		l_index2 := list.New()
		for index := 0; index < N*X; index++ {
		    l_index2 = l_exact[index]
		    e2 := l_index2.Front()
		    d2 := distance(e2.Value.(map[int]int)[0], e2.Value.(map[int]int)[2])
		    a[index] = d2
		    e2 = e2.Next()
	        }

		l_index1 := list.New()
		for index := 0; index < N*X; index++ {
		    l_index1 = l_proj[index]
		    e1 := l_index1.Front()
		    d1 := distance(e1.Value.(map[int]int)[0], e1.Value.(map[int]int)[2])
		    b[index] = d1
		    e1 = e1.Next()
	        }

		dtw := dtw.Dtw{}

		// optionally set your own distance function
		dtw.DistanceFunction = func(x float64, y float64) float64 {
		    difference := x - y
		    return math.Sqrt(difference * difference)
		}
		dtw.ComputeOptimalPathWithWindow(a, b, 5) // 5 = window size
		path := dtw.RetrieveOptimalPath()
                fmt.Println(path)
	}


}
