<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Variational optimiser test on all tet cube/sphere</description>
    <executable>NekMesh</executable>
    <parameters>cube-sphere.msh test.xml:xml:test -v -m varopti:hyperelastic:maxiter=5</parameters>
    <files>
        <file description="Input File">cube-sphere.msh</file>
    </files>
    <metrics>
        <metric type="regex" id="0">
            <regex>Worst at end: (-?\d+(?:\.\d*)?(?:[eE][+\-]?\d+)?)</regex>
            <matches>
                <match>
                    <field id="0" tolerance="5e-2">7.042539e-01</field>
                </match>
            </matches>
        </metric>
    </metrics>
</test>
