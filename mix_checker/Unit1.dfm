object Form1: TForm1
  Left = 242
  Top = 152
  Width = 1391
  Height = 389
  Caption = #1055#1088#1086#1074#1077#1088#1082#1072' '#1085#1072' '#1087#1077#1088#1077#1084#1077#1096#1080#1074#1072#1085#1080#1077
  Color = clBtnFace
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'MS Sans Serif'
  Font.Style = []
  OldCreateOrder = False
  OnCreate = FormCreate
  PixelsPerInch = 96
  TextHeight = 13
  object Image1: TImage
    Left = 8
    Top = 40
    Width = 1089
    Height = 265
  end
  object Image2: TImage
    Left = 1104
    Top = 40
    Width = 265
    Height = 265
  end
  object Label1: TLabel
    Left = 8
    Top = 16
    Width = 164
    Height = 13
    Caption = #1056#1072#1089#1087#1088#1077#1076#1077#1083#1077#1085#1080#1077' '#1074' '#1087#1083#1086#1089#1082#1086#1089#1090#1080' XY:'
  end
  object Label2: TLabel
    Left = 1104
    Top = 16
    Width = 164
    Height = 13
    Caption = #1056#1072#1089#1087#1088#1077#1076#1077#1083#1077#1085#1080#1077' '#1074' '#1087#1083#1086#1089#1082#1086#1089#1090#1080' YZ:'
  end
  object Button1: TButton
    Left = 8
    Top = 312
    Width = 153
    Height = 33
    Caption = #1054#1090#1082#1088#1099#1090#1100' '#1092#1072#1081#1083
    TabOrder = 0
    OnClick = Button1Click
  end
  object OpenDialog1: TOpenDialog
    Left = 168
    Top = 320
  end
end
