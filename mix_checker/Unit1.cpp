//---------------------------------------------------------------------------

#include <vcl.h>
#pragma hdrstop
#include <stdio.h>

#include "Unit1.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"
TForm1 *Form1;
//---------------------------------------------------------------------------
__fastcall TForm1::TForm1(TComponent* Owner)
        : TForm(Owner)
{
}
//---------------------------------------------------------------------------

void __fastcall TForm1::Button1Click(TObject *Sender)
{
int N = 0, k, r;
double A, L, yz_scale, x_scale, x, y, z;

Image1->Canvas->Brush->Color = clWhite;
Image2->Canvas->Brush->Color = clWhite;
Image1->Canvas->Refresh();
Image1->Canvas->Rectangle(0,1,Image1->Width-1, Image1->Height-1);
Image2->Canvas->Refresh();
Image2->Canvas->Rectangle(0,1,Image2->Width-1, Image2->Height-1);

if (OpenDialog1->Execute()) {
    FILE *data_file = fopen(OpenDialog1->FileName.c_str(), "r");
    fscanf(data_file, "%d\n%le\n%le\n", &N, &A, &L);

    Image1->Canvas->Brush->Color = clRed;
    Image2->Canvas->Brush->Color = clRed;
    x_scale = Image1->Width / L;
    yz_scale = Image1->Height / A;
    r = (x_scale + yz_scale) / 2;

    for (int i = 0; i < N; i++) {
        fscanf(data_file, "%d\n%le %le %le\n", &k, &x, &y, &z);
        if (i % 8 == 0) {
            Image1->Canvas->Ellipse(int(x*x_scale-r), Image1->Height - int(y*yz_scale-r), int(x*x_scale+r), Image1->Height - int(y*yz_scale+r) );
            Image2->Canvas->Ellipse(int(y*yz_scale-r), Image2->Height - int(z*yz_scale-r), int(y*yz_scale+r), Image2->Height - int(z*yz_scale+r) );
        }
        fscanf(data_file, "%le %le %le\n", &x, &y, &z);
    }
}
}
//---------------------------------------------------------------------------

void __fastcall TForm1::FormCreate(TObject *Sender)
{
Image1->Canvas->Refresh();
Image1->Canvas->Rectangle(0,1,Image1->Width-1, Image1->Height-1);
Image2->Canvas->Refresh();
Image2->Canvas->Rectangle(0,1,Image2->Width-1, Image2->Height-1);
}
//---------------------------------------------------------------------------
