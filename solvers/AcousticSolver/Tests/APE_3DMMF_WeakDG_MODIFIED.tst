<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>desc P=300</description>
    <executable>AcousticSolver</executable>
    <parameters>APE_3DMMF_WeakDG_MODIFIED.xml</parameters>
    <files>
        <file description="Session File">APE_3DMMF_WeakDG_MODIFIED.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="p" tolerance="1e-12">1.07698e-10</value>
            <value variable="u" tolerance="1e-12">2.29723e-13</value>
            <value variable="v" tolerance="1e-12">2.28589e-13</value>
            <value variable="w" tolerance="1e-12">2.28354e-13</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="p" tolerance="1e-12">1.01636e-10</value>
            <value variable="u" tolerance="1e-5" >2.37098e-13</value>
            <value variable="v" tolerance="1e-5" >2.36746e-13</value>
            <value variable="w" tolerance="1e-5" >2.31423e-13</value>
        </metric>
    </metrics>
</test>
