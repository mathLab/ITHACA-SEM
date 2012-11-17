<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>3D channel flow, Tet elements, par(3), P=8</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>-I GlobalSysSoln=IterativeMultiLevelStaticCond Tet_channel_m8_par.xml</parameters>
    <processes>3</processes>
    <files>
        <file description="Session File">Tet_channel_m8_par.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">0</value>
            <value variable="v" tolerance="1e-12">0</value>
            <value variable="w" tolerance="1e-12">5.26988e-16</value>
	    <value variable="p" tolerance="1e-12">1.36862e-14</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">1.45672e-15</value>
            <value variable="v" tolerance="1e-12">1.4546e-15</value>
            <value variable="w" tolerance="1e-12">3.01703e-14</value>
	    <value variable="p" tolerance="1e-12">4.72511e-13</value>
        </metric>
    </metrics>
</test>
