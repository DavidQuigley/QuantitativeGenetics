#define _CRT_SECURE_NO_DEPRECATE 1
#include <stdio.h>
#include <string>
#include <algorithm>
#include <iostream>
#include <math.h>
#include <fstream>
#include <vector>
#include <set>
using namespace std;
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/random.hpp>
using namespace boost;
namespace alg = boost::algorithm;
#include <wx/wx.h>
#include <wx/grid.h>
#include "carmen/src/DataStructures.h"
#include "carmen/src/graph.h"
#include "carmen/src/Parser.h"
#include "GridClassifierResults.h"


GridClassifierResults::GridClassifierResults(wxWindow* parent, std::string filename, wxWindowID id, wxSize size) : 
wxGrid(parent, id, wxPoint(-1,-1), size){
    
    this->m_parent = parent;

    this->CreateGrid( 0, 0 );
	this->SetRowLabelSize(0);
    this->SetColLabelSize(0);
    this->AppendCols(4);
    this->AppendRows(12);
    //LIGHT_BLUE = (153, 204, 255)
    //LIGHT_GREY = (255, 255, 204)
    this->set_row_color(0, 153, 204, 255);
    this->set_row_color(1, 255, 255, 204);
    this->set_row_color(6, 153, 204, 255);
    this->set_row_color(7, 255, 255, 204);

    CarmenParser cp(filename);
    std::string lblA_train, lblB_train, lblA_test, lblB_test;
    cp.get_value("Class_A_train", 0, lblA_train);  
	cp.get_value("Class_B_train", 0, lblB_train);
    cp.get_value("Class_A_test", 0, lblA_test);  
	cp.get_value("Class_B_test", 0, lblB_test);

    this->SetCellValue(0, 0, wxString::FromAscii( "TRAIN" )  );
    this->SetCellValue(1, 1, wxString::FromAscii( "All" )  );
    this->SetCellValue(1, 2, wxString::FromAscii( lblA_train.c_str() )  );
    this->SetCellValue(1, 3, wxString::FromAscii( lblB_train.c_str() )  );
    this->SetCellValue(2, 0, wxString::FromAscii( "Accuracy" )  );
    this->SetCellValue(3, 0, wxString::FromAscii( "Sensitivity" )  );
    this->SetCellValue(4, 0, wxString::FromAscii( "Specificity" )  );
    this->SetCellValue(5, 0, wxString::FromAscii( "Predictive Value" )  );

    this->SetCellValue(6, 0, wxString::FromAscii( "TEST" )  );
    this->SetCellValue(7, 1, wxString::FromAscii( "All" )  );
    this->SetCellValue(7, 2, wxString::FromAscii( lblA_test.c_str() )  );
    this->SetCellValue(7, 3, wxString::FromAscii( lblB_test.c_str() )  );
    this->SetCellValue(8, 0, wxString::FromAscii( "Accuracy" )  );
    this->SetCellValue(9, 0, wxString::FromAscii( "Sensitivity" )  );
    this->SetCellValue(10, 0, wxString::FromAscii( "Specificity" )  );
    this->SetCellValue(11, 0, wxString::FromAscii( "Predictive Value" )  );

    std::string acc_all_tr, acc_a_tr, acc_b_tr, acc_all_tst, acc_a_tst, acc_b_tst;
    cp.get_value("Acc", 0, acc_all_tr);  
    cp.get_value("Acc", 1, acc_a_tr);  
    cp.get_value("Acc", 2, acc_b_tr);  
    cp.get_value("Acc", 3, acc_all_tst);  
    cp.get_value("Acc", 4, acc_a_tst);  
    cp.get_value("Acc", 5, acc_b_tst);  
    this->SetCellValue(2, 1, wxString::FromAscii( acc_all_tr.c_str() )  );
    this->SetCellValue(2, 2, wxString::FromAscii( acc_a_tr.c_str() )  );
    this->SetCellValue(2, 3, wxString::FromAscii( acc_b_tr.c_str() )  );
    this->SetCellValue(8, 1, wxString::FromAscii( acc_all_tst.c_str() )  );
    this->SetCellValue(8, 2, wxString::FromAscii( acc_a_tst.c_str() )  );
    this->SetCellValue(8, 3, wxString::FromAscii( acc_b_tst.c_str() )  );
    
    std::string sens_tr, sens_tst, spec_tr, spec_tst;
    cp.get_value("Sens_train", 0, sens_tr);
    cp.get_value("Spec_train", 0, spec_tr);
    cp.get_value("Sens_test", 0, sens_tst);
    cp.get_value("Spec_test", 0, spec_tst);
    this->SetCellValue(3, 1, wxString::FromAscii( sens_tr.c_str() )  );
    this->SetCellValue(4, 1, wxString::FromAscii( spec_tr.c_str() )  );
    this->SetCellValue(9, 1, wxString::FromAscii( sens_tst.c_str() )  );
    this->SetCellValue(10, 1, wxString::FromAscii( spec_tst.c_str() )  );

    std::string pv_a_tr, pv_b_tr, pv_a_tst, pv_b_tst;
    cp.get_value("PPV_train", 0, pv_a_tr);
    cp.get_value("NPV_train", 0, pv_b_tr);
    cp.get_value("PPV_test", 0, pv_a_tst);
    cp.get_value("NPV_test", 0, pv_b_tst);
    this->SetCellValue(5, 2, wxString::FromAscii( pv_b_tr.c_str() )  );
    this->SetCellValue(5, 3, wxString::FromAscii( pv_a_tr.c_str() )  );
    this->SetCellValue(11, 2, wxString::FromAscii( pv_b_tst.c_str() )  );
    this->SetCellValue(11, 3, wxString::FromAscii( pv_a_tst.c_str() )  );
    

    this->EnableEditing(false);
    this->Fit();
	this->m_parent->FitInside();
}


void GridClassifierResults::set_row_color(int row, int r, int g, int b){
    for(int i=0; i<4; i++){
        this->SetCellBackgroundColour(row, i, wxColor(r, g, b));
    }
}
