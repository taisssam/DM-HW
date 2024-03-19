#include "Graph.h"

using namespace std;

int main() {

    Graph Europa("europe.txt");

    //Task 2
    cout << "-----------||| TASK 2 |||-----------" << '\n';
    Europa.GetEdgeCount();
    cout << '\n';
    Europa.GetVertexCount();
    cout << '\n';
    Europa.DegreesCount();
    cout << '\n';
    Europa.EccentricitiesCount();
    cout << '\n';
//    Europa.CountCyclomaticNum();
//    cout << '\n';

    //Task 3 - chromatic number
    cout << "-----------||| TASK 3 |||-----------" << '\n';
    Europa.CountChromaticNum();
    cout << '\n';


    //Task 4 - Maximal cliques
    cout << '\n' << "-----------||| TASK 4 |||-----------" << '\n';
    cout << "Cliques: " << '\n';
    Europa.FindMaximalCliques();

    //Task 7 - articulation points
    cout << '\n' << "-----------||| TASK 7 |||-----------" << '\n';
    Europa.Vcomp();
    cout << '\n';

    Europa.Reset();

    //Task 8 - bridges
    cout << '\n' << "-----------||| TASK 8 |||-----------";
    cout << '\n';
    Europa.Ecomp();
    cout << '\n';

    cout << "----------------------------------";

}