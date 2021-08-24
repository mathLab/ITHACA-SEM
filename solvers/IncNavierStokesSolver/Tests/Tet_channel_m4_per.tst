<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>3D channel flow, Tetrahedral elements, P=4, periodic BCs</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>-I GlobalSysSoln=DirectStaticCond Tet_channel_m4_per.xml</parameters>
    <files>
        <file description="Session File">Tet_channel_m4_per.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">2.169e-15</value>
            <value variable="v" tolerance="1e-12">2.23196e-16</value>
            <value variable="w" tolerance="1e-12">4.49841e-16</value>
            <value variable="p" tolerance="1e-8">1.65856e-14</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">2.10942e-15</value>
            <value variable="v" tolerance="1e-12">3.00716e-16</value>
            <value variable="w" tolerance="1e-12">3.21388e-16</value>
            <value variable="p" tolerance="1e-8">1.02141e-14</value>
        </metric>
    </metrics>
</test>
