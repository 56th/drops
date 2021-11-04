#!/usr/bin/env wolframscript
(* ::Package:: *)

path = Directory[] <> "/" <> $ScriptCommandLine[[2]];
img = ImageCrop@Import@path;
path = StringRiffle[Insert[StringSplit[path, "."], "cropped", -2], "."];
Export[path, img]
