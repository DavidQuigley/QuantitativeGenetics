#define _CRT_SECURE_NO_DEPRECATE 1
#include <time.h>
#include <wx/dialog.h>
#include <wx/sizer.h>
#include <wx/stattext.h>
#include <wx/button.h>
#include <wx/textctrl.h>
#include <wx/panel.h>
#include <wx/statbox.h>
#include <wx/wx.h>
#include <wx/process.h>
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/random.hpp>
namespace alg = boost::algorithm;
#include <fstream>
#include <string>
#include <vector>
#include "carmen/src/DataStructures.h"
#include "carmen/src/Attributes.h"
#include "carmen/src/Rawdata.h"
#include "carmen/src/Parser.h"
#include "Investigation.h"
using namespace boost;
using namespace std;
#include "PanelLimits.h"
#include "DlgProgress.h"
#include "PanelDiscretization.h"
#include "DlgDifference.h"

IMPLEMENT_CLASS( DifferenceDialog, wxDialog )

/*
 * Deactivated discretization because the use case (discretizing genes to create sample classes)
 * seemed too complex.
 */
DifferenceDialog::DifferenceDialog( Investigation* investigation, std::string filename ){
	this->investigation = investigation;
	SetExtraStyle(wxWS_EX_BLOCK_EVENTS|wxDIALOG_EX_CONTEXTHELP);
	wxDialog::Create( NULL, wxID_ANY, wxT("Calculate Differential Expression"), wxDefaultPosition, wxDefaultSize, wxCAPTION|wxRESIZE_BORDER|wxSYSTEM_MENU );
	CreateControls(filename);
	GetSizer()->Fit(this); // This fits the dialog to the minimum size dictated by the sizers
	GetSizer()->SetSizeHints(this); // This ensures that the dialog cannot be sized smaller than the minimum size
	Centre();
	redraw();
}

void DifferenceDialog::CreateControls(std::string rs_filename)
{
	char timebuf[80];
	struct tm* newtime;
	time_t long_time;
	time( &long_time );
	newtime = localtime( &long_time );
	strftime(timebuf, 80, "%Y_%m_%d_%H_%M_%S", newtime);

	wxStaticBoxSizer* static_perm_sizer = new wxStaticBoxSizer(wxVERTICAL, this, wxString::FromAscii("Permutation testing") );

	wxStaticText* lbl_file_name = new wxStaticText( this, wxID_ANY, wxT("File Name:") );
	this->txt_file_name = new wxTextCtrl( this, ID_DIFF_TXT_FILE_NAME, wxString::FromAscii(timebuf), wxDefaultPosition, wxSize(200, 22) );
	this->txt_file_name->ChangeValue(wxString::FromAscii(timebuf));
	wxStaticText* lbl_trim = new wxStaticText( this, wxID_ANY, wxT("Expression trim %:"));
	this->txt_trim = new wxTextCtrl( this, ID_DIFF_TRIM, wxT("5"), wxDefaultPosition, wxSize(70, 22) );
	this->txt_trim->SetValue( wxT("0") );

	wxStaticText* lbl_doperms = new wxStaticText( this, wxID_ANY, wxT("Calculate permutations"));
	this->chk_doperms= new wxCheckBox(this, ID_DIFF_CHK_DOPERMS, wxString::FromAscii(""));
	wxStaticText* lbl_max_p = new wxStaticText( this, wxID_ANY, wxT("Max. P value to report:"));
	this->txt_max_p = new wxTextCtrl( this, ID_DIFF_TXT_MAX_P, wxT("1"), wxDefaultPosition, wxSize(70, 22) );
	this->txt_max_p->SetValue( wxT("1") );
	//wxStaticText* lbl_discretization = new wxStaticText( this, wxID_ANY, wxT("Discretization:"));
	//this->discretization = new PanelDiscretization(this, ID_DISCRETIZATION_PANEL, this->investigation);
	wxStaticText* lbl_perms = new wxStaticText( this, wxID_ANY, wxT("Number of permutations:"));
	this->txt_perms = new wxTextCtrl( this, ID_DIFF_TXT_MAX_P, wxT("1000"), wxDefaultPosition, wxSize(70, 22) );
	this->txt_perms->SetValue( wxT("1000") );

	wxStaticBox* static_a = new wxStaticBox(this, wxID_ANY, wxString::FromAscii("Class A Definition"));
    wxStaticBoxSizer* static_class_a_sizer = new wxStaticBoxSizer(static_a, wxHORIZONTAL);
	wxStaticBox* static_b = new wxStaticBox(this, wxID_ANY, wxString::FromAscii("Class B Definition"));
    wxStaticBoxSizer* static_class_b_sizer = new wxStaticBoxSizer(static_b, wxHORIZONTAL);
	this->limit_A = new PanelLimits(this, ID_DIFF_LIMIT_A, this->investigation);
	this->limit_B = new PanelLimits(this, ID_DIFF_LIMIT_B, this->investigation);
	this->limit_A->disable_gene_limits();
	this->limit_B->disable_gene_limits();
	if(investigation->is_mac){
		wxFont font(wxNORMAL_FONT->GetPointSize(), wxFONTFAMILY_SWISS, wxNORMAL, wxFONTWEIGHT_BOLD);
		static_a->SetFont( font );
		static_b->SetFont( font );
		static_perm_sizer->GetStaticBox()->SetFont( font );
	}
	if( rs_filename.size()>0 ){
		CarmenParser cp(rs_filename);
		std::string class_a, class_b;
		cp.get_value(std::string("Class_A"), 0, class_a);
		cp.get_value(std::string("Class_B"), 0, class_b);
		this->limit_A->set_limits( class_a );
		this->limit_B->set_limits( class_b );
	}
	wxSizer* sizer_ok = CreateButtonSizer(wxOK|wxCANCEL);
	Connect(wxID_OK, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(DifferenceDialog::OnClickOk));
	Connect(ID_DIFF_LIMIT_A, wxEVT_COMMAND_TEXT_UPDATED, wxCommandEventHandler(DifferenceDialog::OnUpdateLabel));
	Connect(ID_DIFF_LIMIT_B, wxEVT_COMMAND_TEXT_UPDATED, wxCommandEventHandler(DifferenceDialog::OnUpdateLabel));
	Connect(ID_DIFF_CHK_DOPERMS, wxEVT_COMMAND_CHECKBOX_CLICKED, wxCommandEventHandler(DifferenceDialog::OnCheckDoPerms));

	wxBoxSizer* sizer_top = new wxBoxSizer(wxVERTICAL);

	wxFlexGridSizer* sizer_upper = new wxFlexGridSizer(2, 2, 0, 0);
	sizer_upper->Add(lbl_file_name, 0, FLAGS, BORDER_PXL);
	sizer_upper->Add(this->txt_file_name, 0, FLAGS, BORDER_PXL);
	sizer_upper->Add(lbl_trim, 0, FLAGS, BORDER_PXL);
	sizer_upper->Add(this->txt_trim, 0, FLAGS, BORDER_PXL);
	sizer_top->Add(sizer_upper, 0, FLAGS, BORDER_PXL);


	wxFlexGridSizer* sizer_perm = new wxFlexGridSizer(3, 2, 0, 0);
	sizer_perm->Add(lbl_doperms, 0, FLAGS, BORDER_PXL);
	sizer_perm->Add(this->chk_doperms, 0, FLAGS, BORDER_PXL);
	sizer_perm->Add(lbl_max_p, 0, FLAGS, BORDER_PXL);
	sizer_perm->Add(this->txt_max_p, 0, FLAGS, BORDER_PXL);
	sizer_perm->Add(lbl_perms, 0, FLAGS, BORDER_PXL);
	sizer_perm->Add(this->txt_perms, 0, FLAGS, BORDER_PXL);
	static_perm_sizer->Add(sizer_perm, 0, FLAGS, BORDER_PXL);
	sizer_top->Add(static_perm_sizer, 0, FLAGS, BORDER_PXL);
	//sizer_upper->Add(lbl_discretization, 0, FLAGS, BORDER_PXL);
	//sizer_upper->Add(this->discretization, 0, FLAGS, BORDER_PXL);

	wxBoxSizer* sizer_limits = new wxBoxSizer(wxHORIZONTAL);
	static_class_a_sizer->Add(this->limit_A, 0, FLAGS, BORDER_PXL);
	static_class_b_sizer->Add(this->limit_B, 0, FLAGS, BORDER_PXL);
	sizer_limits->Add(static_class_a_sizer, 0, FLAGS, BORDER_PXL);
	sizer_limits->Add(static_class_b_sizer, 0, FLAGS, BORDER_PXL);
	sizer_top->Add(sizer_limits, 0, FLAGS_RIGHT, BORDER_PXL);

	sizer_top->Add(sizer_ok, 0, FLAGS_RIGHT, BORDER_PXL);

	this->SetSizer(sizer_top);
	redraw();
}


void DifferenceDialog::OnCheckDoPerms( wxCommandEvent& evt ){
	redraw();
}

void DifferenceDialog::OnUpdateLabel( wxCommandEvent& evt){
	redraw();
}


void DifferenceDialog::OnClickOk( wxCommandEvent& WXUNUSED(event) ){
	std::string d( this->investigation->fn_expr );
	std::string f( this->investigation->fn_sa );
	std::string g( this->investigation->fn_ga );
	std::string a = this->limit_A->get_limits();
	std::string b = this->limit_B->get_limits();
	std::string p( this->txt_max_p->GetValue().ToAscii() );
    std::string n( this->txt_perms->GetValue().ToAscii() );
	std::string o( this->txt_file_name->GetValue().ToAscii() );
	std::string t( this->txt_trim->GetValue().ToAscii() );
	if( a.size()==0 || b.size()==0){
		wxMessageBox(wxString::FromAscii("Please choose limits for class A and class B"), _T("Not Ready"), wxOK | wxICON_EXCLAMATION, this);
		return;
	}
	if( o.size()==0 ){
		wxMessageBox(wxString::FromAscii("Please enter a name for the difference file that will be generated"), _T("Not Ready"), wxOK | wxICON_EXCLAMATION, this);
		return;
	}
	if( ! alg::find_first( o, ".difference") )
		o = o + ".difference";
	std::stringstream ss;

	if(this->investigation->is_mac)
		ss << this->investigation->dir_results << '/' << o;
	else
		ss << this->investigation->dir_results << '\\' << o;
	std::string attempted_filename( ss.str() );
	float p_f;
	int n_i;
    if( this->chk_doperms->IsChecked() ){
		try{ p_f = boost::lexical_cast<float>(p); }
		catch( ... ){
			wxMessageBox(wxString::FromAscii("p-value limit must be a real-valued number between 0 and 1"), _T("Not Ready"), wxOK | wxICON_EXCLAMATION, this);
			return;
		}
		if( p_f < 0 || p_f > 1 ){
			wxMessageBox(wxString::FromAscii("p-value limit must be a real-valued number between 0 and 1"), _T("Not Ready"), wxOK | wxICON_EXCLAMATION, this);
			return;
		}
		try{ n_i = boost::lexical_cast<int>(n); }
	    catch( ... ){
			wxMessageBox(wxString::FromAscii("Number of permutations must be a positive whole number"), _T("Not Ready"), wxOK | wxICON_EXCLAMATION, this);
			return;
		}
	}
	float t_f;
	try{
		t_f = boost::lexical_cast<float>(t);
	}
	catch( ... ){
		wxMessageBox(wxString::FromAscii("Trim percentage must be a real-valued number between 0 and 100"), _T("Not Ready"), wxOK | wxICON_EXCLAMATION, this);
		return;
	}
	if( t_f < 0 || t_f > 100 ){
		wxMessageBox(wxString::FromAscii("Trim percentage must be a real-valued number between 0 and 100"), _T("Not Ready"), wxOK | wxICON_EXCLAMATION, this);
		return;
	}

	std::string y( this->investigation->ga->get_gene_name_column() );

	std::vector<std::string> cmd;
	cmd.push_back( "difference" );
	cmd.push_back("-d" + d );
	cmd.push_back("-f" + f );
	cmd.push_back("-g" + g );
    cmd.push_back("-a" + a );
    cmd.push_back("-b" + b );
    if( this->chk_doperms->IsChecked() ){
    	cmd.push_back("-p" + p );
    	cmd.push_back("-n" + n );
    }
    else{
    	cmd.push_back("-p1" );
    	cmd.push_back("-n1" );
    }
    cmd.push_back("-t" + t );
    cmd.push_back("-y" + y );
	cmd.push_back("-vT");
    cmd.push_back("-o" + attempted_filename );
	//int d_idx = this->discretization->get_discretization_code();
    //if( d_idx==DISC_NONE )
	//	cmd.push_back("-mnone");
	//else if( d_idx==DISC_SD ){
	//	cmd.push_back("-mSD");
	//}
    //else if( d_idx==DISC_MAD ){
	//	cmd.push_back("-mMAD");
	//}
    //else if( d_idx==DISC_ABS ){
	//	cmd.push_back("-mabs");
	//}
    //else if( d_idx==DISC_PER ){
	//	cmd.push_back("-mper");
	//}
	//if( d_idx != DISC_NONE ){
	//	std::string lower, upper;
	//	this->discretization->get_discretization_bounds(lower, upper);
	//	cmd.push_back("-l" + lower );
    //    cmd.push_back("-u" + upper );
	//}
	ProgressDialog dialog(this, cmd, this->investigation, true);
	if (dialog.ShowModal() == wxID_OK){
        this->new_filename = attempted_filename;
		this->AcceptAndClose();
	}
}


void DifferenceDialog::redraw(){
	if( this->chk_doperms->IsChecked() ){
		this->txt_perms->Enable(true);
		this->txt_max_p->Enable(true);
        this->txt_perms->SetBackgroundColour(*wxWHITE);
        this->txt_max_p->SetBackgroundColour(*wxWHITE);
	}
	else{
		this->txt_perms->Enable(false);
		this->txt_max_p->Enable(false);
        this->txt_perms->SetBackgroundColour(*wxLIGHT_GREY);
        this->txt_max_p->SetBackgroundColour(*wxLIGHT_GREY);
	}
	//std::string lim_A = this->limit_A->get_limits();
	//std::string lim_B = this->limit_A->get_limits();
	//if( boost::algorithm::find_first(lim_A, "gene:") || boost::algorithm::find_first(lim_B, "gene:") )
	//	this->discretization->Enable(true);
	//else
	//	this->discretization->Enable(false);
}
