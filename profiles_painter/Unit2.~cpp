//---------------------------------------------------------------------------

#include <vcl.h>
#pragma hdrstop
#include <stdio.h>
#include <math.h>
#include "Unit1.h"
#include "Unit2.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"
TForm2 *Form2;
//---------------------------------------------------------------------------
__fastcall TForm2::TForm2(TComponent* Owner)
        : TForm(Owner)
{
}
//---------------------------------------------------------------------------

void __fastcall TForm2::Button1Click(TObject *Sender)
{
        if (ColorDialog1->Execute()) {
                Shape1->Brush->Color = ColorDialog1->Color;
        }
}
//---------------------------------------------------------------------------

void __fastcall TForm2::Button2Click(TObject *Sender)
{
       if (OpenDialog1->Execute()) {
                Edit1->Text = OpenDialog1->FileName;
        }
}
//---------------------------------------------------------------------------
void __fastcall TForm2::Button3Click(TObject *Sender)
{
        double IMG_Y;
        double IMG_Z;

        int start_x = 50, start_y = 50;

        FILE *file = fopen(Edit1->Text.c_str(), "r");

        while (fscanf(file, "%le\n", &IMG_Y) + fscanf(file, "%le\n", &IMG_Z) > 0) {
                int y = ceill(IMG_Y*20.0L);
                int z = ceill(IMG_Z*20.0L);
                Form1->Image1->Canvas->Pixels[ start_x + y ][ Form1->Image1->Height - start_y - z ] = Shape1->Brush->Color;
        }

        Form2->Hide();
}
//---------------------------------------------------------------------------

void __fastcall TForm2::Shape1MouseDown(TObject *Sender,
      TMouseButton Button, TShiftState Shift, int X, int Y)
{
        if (ColorDialog1->Execute()) {
                Shape1->Brush->Color = ColorDialog1->Color;
        }        
}
//---------------------------------------------------------------------------
