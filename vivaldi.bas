'#Language "WWB-COM"

Option Explicit

Sub Main


Dim I As Long, myString
Dim iMyFreeFile As Integer

 iMyFreeFile = FreeFile
 'reset any file accidentally left open


OpenFile("D:\Emad Arasteh\thesis\summer 2019\opt\PSO_optimization\TR2-.cst")
 ' Open file for input.
SendKeys "{F7}"

Open "D:\Emad Arasteh\thesis\summer 2019\opt\PSO_optimization\TR2-\Export\faz1.txt" For Input As #1


 Input #1, myString
storeparameter("phi1",myString)
Close #1
SendKeys "{F7}"


Open "D:\Emad Arasteh\thesis\summer 2019\opt\PSO_optimization\TR2-\Export\faz2.txt" For Input As #2

 Input #2, myString
storeparameter("phi2",myString)
Close #2
SendKeys "{F7}"

Open "D:\Emad Arasteh\thesis\summer 2019\opt\PSO_optimization\TR2-\Export\faz3.txt" For Input As #3

 Input #3, myString
storeparameter("phi3",myString)
Close #3
SendKeys "{F7}"


Open "D:\Emad Arasteh\thesis\summer 2019\opt\PSO_optimization\TR2-\Export\faz4.txt" For Input As #4

 Input #4, myString
storeparameter("phi4",myString)
Close #4
SendKeys "{F7}"

Open "D:\Emad Arasteh\thesis\summer 2019\opt\PSO_optimization\TR2-\Export\faz5.txt" For Input As #5

 Input #5, myString
storeparameter("phi5",myString)
Close #5
SendKeys "{F7}"


Open "D:\Emad Arasteh\thesis\summer 2019\opt\PSO_optimization\TR2-\Export\faz6.txt" For Input As #6

 Input #6, myString
storeparameter("phi6",myString)
Close #6
SendKeys "{F7}"

Open "D:\Emad Arasteh\thesis\summer 2019\opt\PSO_optimization\TR2-\Export\faz7.txt" For Input As #7

 Input #7, myString
storeparameter("phi7",myString)
Close #7
SendKeys "{F7}"


Open "D:\Emad Arasteh\thesis\summer 2019\opt\PSO_optimization\TR2-\Export\faz8.txt" For Input As #8

 Input #8, myString
storeparameter("phi8",myString)
Close #8
SendKeys "{F7}"








For I = 1 To 8
   SendKeys "{F7}"
Next I
'Reset


SendKeys "{F7}"






Solver.Start
SendKeys "{F7}"
save

SendKeys "{F7}"







End Sub
