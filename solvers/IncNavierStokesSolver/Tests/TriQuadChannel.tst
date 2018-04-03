<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Tri and Quad element channel Flow 2D</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>TriQuadChannel.xml</parameters>
    <files>
        <file description="Session File">TriQuadChannel.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">5.37124e-14</value>
            <value variable="v" tolerance="1e-12">2.48979e-14</value>
	    <value variable="p" tolerance="5e-12">0.0</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">1.47216e-13</value>
            <value variable="v" tolerance="1e-12">2.96768e-14</value>
	    <value variable="p" tolerance="1e-12">7.59719e-13</value>
        </metric>
    </metrics>
</test>


