<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Post-process Semtex flow field of Kovasznay flow</description>
    <executable python="true">kovas2.py</executable>
    <parameters></parameters>
    <files>
        <file description="Session File">kovas2.xml</file>
        <file description="Field file">kovas2.sem.fld</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="x" tolerance="1e-12">1.0854</value>
            <value variable="y" tolerance="1e-12">0.626657</value>
            <value variable="z" tolerance="1e-12">6.37932</value>
            <value variable="u" tolerance="1e-12">3.65663</value>
            <value variable="v" tolerance="1e-12">0.313122</value>
            <value variable="w" tolerance="1e-12">0</value>
            <value variable="p" tolerance="1e-8">1.06858</value>
        </metric>
    </metrics>
</test>
