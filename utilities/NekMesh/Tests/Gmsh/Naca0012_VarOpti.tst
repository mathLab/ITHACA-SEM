<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Variational optimiser test on a 2D all-triangle NACA0012 case</description>
    <executable>NekMesh</executable>
    <parameters>Naca0012_VarOpti.msh test.xml:xml:test -v -m varopti:linearelastic:numthreads=2:maxiter=5</parameters>
    <files>
        <file description="Input File">Naca0012_VarOpti.msh</file>
    </files>
    <metrics>
        <metric type="regex" id="0">
            <regex>Worst at end: (-?\d+(?:\.\d*)?(?:[eE][+\-]?\d+)?)</regex>
            <matches>
                <match>
                    <field id="0">8.725848e-01</field>
                </match>
            </matches>
        </metric>
    </metrics>
</test>
