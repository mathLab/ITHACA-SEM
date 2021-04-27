<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> Fitst part of wallNormalData module (projection a point to the boudndary) </description>
    <executable>FieldConvert</executable>
    <parameters> -f -e -m wallNormalData:bnd=0:xorig="0,0":projDir="1,1":nptsH=3:distH=0.002 wallNormalData.xml wallNormalData.fld wallNormalData.pts</parameters>
    <files>
        <file description="Session File">wallNormalData.xml</file>
        <file description="Session File">wallNormalData.fld</file>
    </files>
    <metrics>
        <metric type="Linf" id="1">
            <value variable="x" tolerance="1e-6">0.354968</value>
            <value variable="y" tolerance="1e-6">0.354968</value>
        </metric>
    </metrics>
</test>
