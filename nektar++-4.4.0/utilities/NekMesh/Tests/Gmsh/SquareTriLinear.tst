<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Gmsh linear tri square, order 9 output</description>
    <executable>NekMesh</executable>
    <parameters>-m jac:list SquareTriLinear.msh SquareTri.xml:xml:test:order=9</parameters>
    <files>
        <file description="Input File">SquareTriLinear.msh</file>
    </files>
    <metrics>
        <metric type="regex" id="1">
            <regex>^Total negative Jacobians: (\d+)</regex>
            <matches>
                <match>
                    <field id="0">0</field>
                </match>
            </matches>
        </metric>
    </metrics>
</test>
