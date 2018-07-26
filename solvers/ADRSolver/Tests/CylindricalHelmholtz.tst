<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> 2D Helmholtz with cylindrical periodicity P=5 </description>
    <executable>ADRSolver</executable>
    <parameters>CylindricalHelmholtz.xml</parameters>
    <files>
        <file description="Session File">CylindricalHelmholtz.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-10"> 4.95977e-06 </value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-10"> 0.000144863 </value>
        </metric>
    </metrics>
</test>
