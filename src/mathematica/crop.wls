#!/usr/bin/env wolframscript
(* ::Package:: *)

path = Directory[] <> "/" <> $ScriptCommandLine[[2]];
If[Length[$ScriptCommandLine] == 3, s = ToExpression[$ScriptCommandLine[[3]]]; img = ImageCrop[Import@path, {s, s}], img = ImageCrop@Import@path];
path = StringRiffle[Insert[StringSplit[path, "."], "cropped", -2], "."];
Export[path, img]
