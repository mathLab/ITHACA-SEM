<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Helmholtz 3D CG for Tet iterative SL static global par(3)</description>
    <executable>Helmholtz3D</executable>
    <parameters>--use-scotch Helmholtz3D_Tet_iter_global.xml</parameters>
    <processes>3</processes>
    <files>
        <file description="Session File">Helmholtz3D_Tet_iter_global.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-8">5.48966e-05</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-8">0.000344037</value>
        </metric>
    </metrics>
</test>


