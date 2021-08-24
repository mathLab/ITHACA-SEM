<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Add composite ID to field</description>
    <executable>FieldConvert</executable>
    <parameters> -f -e -m addcompositeid compositeid.xml compositeid.fld</parameters>
    <files>
        <file description="Session File">compositeid.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="x" tolerance="1e-6">3</value>
            <value variable="y" tolerance="1e-6">1</value>
            <value variable="compositeID" tolerance="1e-12">2.23607</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="x" tolerance="1e-6">3</value>
            <value variable="y" tolerance="1e-6">1</value>
            <value variable="compositeID" tolerance="1e-12">2</value>
        </metric>
    </metrics>
</test>
