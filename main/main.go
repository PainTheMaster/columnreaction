package main

import "PainTheMaster/columnreaction/column"

func main() {
	var a, b column.AnalyteInCol

	totalT := 15.0 * 60.0                 // s
	divTime := 100000                     //
	timeFrac := totalT / float64(divTime) //s

	a.Init(1.0, //M iniConc float64
		10.0/100.0,       //m colLen float64,
		1000,             //division int,
		10.0e-2/(1.0*60), // m/s flowRate float64
		0.1,              //relRate float64
		4.0e-6,           //diffRate float64,
		5.0e-4,           // decompRate float64
		totalT,           // float64
		timeFrac,         //s float64
	)

	b.Init(0.2, //M iniConc float64
		10.0/100.0,       //m colLen float64,
		1000,             //division int,
		10.0e-2/(1.0*60), // m/s flowRate float64
		0.15,             //relRate float64
		4.0e-6,           //diffRate float64,
		0.0,              // decompRate float64
		totalT,           // float64
		timeFrac,         //s float64
	)

	for i := 1; i <= divTime; i++ {
		a.Diffuse()
		b.Diffuse()
		a.React(&b)
		a.Transl(i)
		b.Transl(i)
	}

	a.Output("a.csv")
	b.Output("b.csv")

}
