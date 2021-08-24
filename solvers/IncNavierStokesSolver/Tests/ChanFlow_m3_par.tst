<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>3D channel flow, 2D, par(2), P=8</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>--use-scotch ChanFlow_m3_par.xml</parameters>
    <processes>2</processes>
    <files>
        <file description="Session File">ChanFlow_m3_par.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-6">8.42518e-14</value>
            <value variable="v" tolerance="1e-6">4.35393e-13</value>
	    <value variable="p" tolerance="1e-6">6.56835e-11</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-6">6.34659e-13</value>
            <value variable="v" tolerance="1e-6">6.83331e-13</value>
	    <value variable="p" tolerance="1e-6">1.68574e-10</value>
        </metric>
    </metrics>
</test>


