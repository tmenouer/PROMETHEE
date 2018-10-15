package main

import (
       "math"
       "math/rand"
       "fmt"
       "sync"
       "time"
       //"container/list"
       "strconv"
)

/*
 * LSH ideas from https://github.com/ekzhu/lsh
 */

// Value is an index into the input dataset.
type hashTableBucket []string

// Point is a vector in the L2 metric space.
type Point map[int]int
type PointP []float64

// Dot returns the dot product of two points.
func (p Point) Dot(q Point) float64 {
	s := 0
	for i := 0; i < len(p); i = i + 2 {
		s += p[i] * q[i]
	}
	return float64(s)
}

// L2 returns the L2 distance of two points.
func (p Point) L2(q Point) float64 {
	s := 0
	for i := 0; i < len(p); i = i + 2 {
		d := p[i] - q[i]
		s += d * d
	}
	return math.Sqrt(float64(s))
}

type lshParams struct {
	// Dimensionality of the input data.
	dim int
	// Number of hash tables.
	l int
	// Number of hash functions for each table.
	m int
	// Shared constant for each table.
	w float64

	// Hash function params for each (l, m).
	a [][]Point
	b [][]float64
}

// NewLshParams initializes the LSH settings.
func newLshParams(dim, l, m int, w float64) *lshParams {
	// Initialize hash params.
	a := make([][]Point, l)
	b := make([][]float64, l)
	random := rand.New(rand.NewSource(time.Now().UTC().UnixNano()))
	for i := range a {
		a[i] = make([]Point, m)
		b[i] = make([]float64, m)
		for j := range a[i] {
			a[i][j] = make(Point, dim)
			for d := 0; d < dim; d++ {
			        // Attention au cast. Modification de Christophe Cerin
				a[i][j][d] = int(random.NormFloat64())
			}
			b[i][j] = random.Float64() * float64(w)
		}
	}
	return &lshParams{
		dim: dim,
		l:   l,
		m:   m,
		a:   a,
		b:   b,
		w:   w,
	}
}

// Key is a way to index into a table.
type hashTableKey []int

// Hash returns all combined hash values for all hash tables.
func (lsh *lshParams) hash(point Point) []hashTableKey {
	hvs := make([]hashTableKey, lsh.l)
	for i := range hvs {
		s := make(hashTableKey, lsh.m)
		for j := 0; j < lsh.m; j++ {
			hv := (point.Dot(lsh.a[i][j]) + lsh.b[i][j]) / lsh.w
			s[j] = int(math.Floor(hv))
		}
		hvs[i] = s
	}
	return hvs
}

type basicHashTableKey string

type hashTable map[basicHashTableKey]hashTableBucket

// BasicLsh implements the original LSH algorithm for L2 distance.
type BasicLsh struct {
	*lshParams
	// Hash tables.
	tables []hashTable
}

// NewBasicLsh creates a basic LSH for L2 distance.
// dim is the diminsionality of the data, l is the number of hash
// tables to use, m is the number of hash values to concatenate to
// form the key to the hash tables, w is the slot size for the
// family of LSH functions.
func NewBasicLsh(dim, l, m int, w float64) *BasicLsh {
	tables := make([]hashTable, l)
	for i := range tables {
		tables[i] = make(hashTable)
	}
	return &BasicLsh{
		lshParams: newLshParams(dim, l, m, w),
		tables:    tables,
	}
}

func (index *BasicLsh) toBasicHashTableKeys(keys []hashTableKey) []basicHashTableKey {
	basicKeys := make([]basicHashTableKey, index.l)
	for i, key := range keys {
		s := ""
		for _, hashVal := range key {
			s += fmt.Sprintf("%.16x", hashVal)
		}
		basicKeys[i] = basicHashTableKey(s)
	}
	return basicKeys
}

// Insert adds a new data point to the LSH.
// id is the unique identifier for the data point.
func (index *BasicLsh) Insert(point Point, id string) {
	// Apply hash functions
	hvs := index.toBasicHashTableKeys(index.hash(point))
	//index.toBasicHashTableKeys(index.hash(point))
	//index.hash(point)
	// Insert key into all hash tables
	/**/
	var wg sync.WaitGroup
	wg.Add(len(index.tables))
	/**/
	for i := range index.tables {
		hv := hvs[i]
		table := index.tables[i]
		/* */
		go
		func(table hashTable, hv basicHashTableKey) {
		/**/
		if _, exist := table[hv]; !exist {
				table[hv] = make(hashTableBucket, 0)
		}
		table[hv] = append(table[hv], id)
		/**/
			wg.Done()
		}(table, hv)
		/**/
	}
	/**/
	wg.Wait()
	/**/
}

// Insert adds a new data point to the LSH.
// id is the unique identifier for the data point.
func (index *BasicLsh) Insert_Seq(point Point, id string) {
	// Apply hash functions
	hvs := index.toBasicHashTableKeys(index.hash(point))
	//index.toBasicHashTableKeys(index.hash(point))
	//index.hash(point)
	// Insert key into all hash tables
	/**/
	//var wg sync.WaitGroup
	//wg.Add(len(index.tables))
	/**/
	//fmt.Println("**********")
	for i := range index.tables {
		hv := hvs[i]
		table := index.tables[i]
		/* */
		//go
		//func(table hashTable, hv basicHashTableKey) {
		/**/
		//fmt.Println("*********")
		if _, exist := table[hv]; !exist {
				table[hv] = make(hashTableBucket, 0)
		}
		//fmt.Println("*********")
		table[hv] = append(table[hv], id)
		/**/
		//	wg.Done()
		//}(table, hv)
		/**/
	}
	//fmt.Println("------------------")
	/**/
	//wg.Wait()
	/**/
}

// Query finds the ids of approximate nearest neighbour candidates,
// in un-sorted order, given the query point,
func (index *BasicLsh) Query(q Point) []string {
	// Apply hash functions
	hvs := index.toBasicHashTableKeys(index.hash(q))
	// Keep track of keys seen
	seen := make(map[string]bool)
	for i, table := range index.tables {
		if candidates, exist := table[hvs[i]]; exist {
			for _, id := range candidates {
				if _, exist := seen[id]; exist {
					continue
				}
				seen[id] = true
			}
		}
	}
	// Collect results
	ids := make([]string, 0, len(seen))
	for id := range seen {
		ids = append(ids, id)
	}
	return ids
}

// randomPoints returns a slice of point vectors,
// each element of every point vector is drawn from a uniform
// distribution over [0, max)
func randomPoints(n, Dim1, Dim2  int) []Point {
     MyTab := make([]Point, n) // Point == Point 
     for i := range MyTab {
		MyTab[i] = make(Point)
		MyTab[i][0] = rand.Intn(Dim1) //Value in x for function 1
		MyTab[i][1] = 0               // O if minimization - 1 if maximization
		MyTab[i][2] = rand.Intn(Dim2) //value in y for function 2
		MyTab[i][3] = 0               // O if minimization - 1 if maximization
		MyTab[i][4] = i               //Id value
	}
	return MyTab
	
}

func main() {
	//
	// Parameters of the experiment
	//
	powof2 := int(math.Pow(2.0, 20.0)) // Max size ;
	MinSize := 4096                    // Min size = 4096 by default
	Dim1 := 100                   // random values between 0 and Dim1
	Dim2 := 100                   // random values between 0 and Dim2
	const N int = 15                   // Number of runs
	const NN int = int(N) + 1          // Number of runs + 1

	fmt.Println("Begin...")
	rand.Seed(time.Now().UTC().UnixNano())
	var result_3way [NN]map[int]float64
	var result_kung [NN]map[int]float64

	for j := 0; j < int(N)+1; j++ {
		result_3way[j] = make(map[int]float64)
		result_kung[j] = make(map[int]float64)
	}

	for j := 0; j < int(N); j++ {
		for taille := MinSize; taille <= powof2; taille = taille * 2 {

			//l := list.New()

			fmt.Println("Problem size ", taille)
			
			fmt.Println("Building of a new LSH dict")
			// The m (Log factor) parameter is too high to get performance
			//lsh := NewBasicLsh(2, 1, int(math.Log2(float64(taille))), 5.0)
			// m = 2 provides with good results
			lsh := NewBasicLsh(2, 1, int(2), 5.0)
			fmt.Println("Build a random sample and make buckets")
			
			points := randomPoints(taille, Dim1, Dim2)
			
			start := time.Now()
			insertedKeys := make([]string, taille)
			//start := time.Now()
			for i, p := range points {
			    //fmt.Println(i," ",p)
			    lsh.Insert_Seq(p, strconv.Itoa(i))
			    insertedKeys[i] = strconv.Itoa(i)
			}
			elapsed := time.Since(start)
			fmt.Printf("Insertion of points in buckets takes %s\n", elapsed)
			
			// Use the inserted points as queries, and
			// verify that we can get back each query itself
			/*
			for i, key := range insertedKeys {
			    found := false
			    for _, foundKey := range lsh.Query(points[i]) {
			    	if foundKey == key {
				   found = true
				}
			    }
			    if !found {
			       fmt.Println("Query fail")
			    }
			}
			*/

		}
	}


	// We print a summary of results
	/*
	fmt.Println("-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-")

	fmt.Println("Mean execution time (N experiments) for LSH")
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

	points := [][]float64{{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}}
	points_orig := [][]float64{{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}}
	*/
	/*
		points := [][]float64{{0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}}
		points_orig := [][]float64{{0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}}
	*/
	/*
	i := 0
	max1 := 0.0
	for key := MinSize; key <= powof2; key *= 2 {
		points[0][i] = float64(key)
		points[1][i] = result_3way[N][key] / float64(N)
		i++
		if result_3way[N][key]/float64(N) > max1 {
			max1 = result_3way[N][key] / float64(N)
		}
	}

	i = 0
	max2 := 0.0
	for key := MinSize; key <= powof2; key *= 2 {
		points_orig[0][i] = float64(key)
		points_orig[1][i] = result_kung[N][key] / float64(N)
		i++
		if result_kung[N][key]/float64(N) > max2 {
			max2 = result_kung[N][key] / float64(N)
		}
	}

	var maxmax = 0.0
	if max2 > max1 {
		maxmax = max2
	} else {
		maxmax = max1
	}
	*/
	//fmt.Println(maxmax)
	//fmt.Println(points_orig)
	/*
	fmt.Println("Print on the screen a Gnuplot program to execute on the command line")
	str := "\"front" + strconv.FormatFloat(float64(time.Now().UTC().UnixNano()), 'f', -1, 64) + ".svg\""
	fmt.Println("Picture name: ", str)
	fmt.Println("set terminal svg linewidth 1 font 'Arial,30' size 1280,960 ;")
	fmt.Println("set logscale x 2 ;set format x \"2^{%L}\" ; set xtics center logscale;")
	fmt.Println("set output ", str)
	fmt.Printf("set yrange [0:%d]\n", int(maxmax*1.115))
	fmt.Println("set xlabel \"Input size (tuples)\"")
	fmt.Println("set ylabel \"Time (ms)\"")
	fmt.Println("plot '-' using 1:2 with linespoints title \"3way\", '-' using 1:2 with linespoints title \"Kung\" ")
	for key := MinSize; key <= powof2; key *= 2 {
		fmt.Println(key, " ", result_3way[N][key]/float64(N))
	}
	fmt.Println("e")
	for key := MinSize; key <= powof2; key *= 2 {
		fmt.Println(key, " ", result_kung[N][key]/float64(N))
	}
	fmt.Println("e")
	*/
}
