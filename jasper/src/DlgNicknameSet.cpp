#include <string>
#include <vector>
#include <set>
#include <wx/dialog.h>
#include <wx/sizer.h>
#include <wx/stattext.h>
#include <wx/choice.h>
#include <wx/button.h>
#include <wx/arrstr.h>
#include <wx/textctrl.h>
#include <wx/wx.h>
#include <boost/algorithm/string.hpp>
namespace alg = boost::algorithm;
#include <boost/tokenizer.hpp>
#include <boost/random.hpp>
using namespace boost;
using namespace std;
#include "DlgNicknameSet.h"

IMPLEMENT_CLASS( NicknameSetDialog, wxDialog )

NicknameSetDialog::NicknameSetDialog( ){
}

NicknameSetDialog::NicknameSetDialog( wxWindow* parent, wxWindowID id, const wxString& caption, const wxPoint& pos, const wxSize& size, long style )
{
	Create(parent, id, caption, pos, size, style);
}

bool NicknameSetDialog::Create( wxWindow* parent, wxWindowID id, const wxString& caption, const wxPoint& pos, const wxSize& size, long style )
{
	SetExtraStyle(wxWS_EX_BLOCK_EVENTS|wxDIALOG_EX_CONTEXTHELP);
	if (!wxDialog::Create( parent, id, caption, pos, size, style ))
		return false;
	CreateControls();
	GetSizer()->Fit(this); // This fits the dialog to the minimum size dictated by the sizers
	GetSizer()->SetSizeHints(this); // This ensures that the dialog cannot be sized smaller than the minimum size
	Centre(); 
	return true;
}


void NicknameSetDialog::CreateControls()
{
	// CONTROLS
	wxStaticText* lbl_nickname = new wxStaticText( this, wxID_ANY, wxT("Nickname:") );
	wxTextCtrl* txt_nickname = new wxTextCtrl( this, ID_TXT_NICKNAME );
	
	wxButton* btn_ok = new wxButton( this, wxID_OK, wxT("&Ok"), wxDefaultPosition, wxDefaultSize, 0 );
	wxButton* btn_cancel = new wxButton( this, wxID_CANCEL, wxT("&Cancel"), wxDefaultPosition, wxDefaultSize, 0 );
	btn_ok->SetDefault();

	Connect(wxID_OK, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(NicknameSetDialog::OnClickOk));

	// POSITIONING
	wxFlexGridSizer* sizer_upper = new wxFlexGridSizer(2, 2, 0, 0);
	sizer_upper->Add(lbl_nickname, 0, FLAGS, BORDER_PXL);
	sizer_upper->Add(txt_nickname, 0, FLAGS, BORDER_PXL);

	sizer_upper->Add(btn_ok, 0, FLAGS, BORDER_PXL);
	sizer_upper->Add(btn_cancel, 0, FLAGS, BORDER_PXL);

	wxBoxSizer* sizer_top = new wxBoxSizer(wxVERTICAL);
	sizer_top->Add(sizer_upper, 0, FLAGS, BORDER_PXL);
	this->SetSizer(sizer_top);
}


void NicknameSetDialog::OnClickOk( wxCommandEvent& event ){
	wxString newnick = ((wxTextCtrl*)FindWindow(ID_TXT_NICKNAME))->GetValue();
	this->new_nickname = std::string( newnick.ToAscii() );
	if( this->new_nickname.size() > 0 )
		this->AcceptAndClose();
	else
		wxMessageBox(wxString::FromAscii("ERROR: Nickname must not be blank"), _T("Error"), wxOK | wxICON_INFORMATION, this);
}
