<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Gmsh (v4.0) high-order tet cube, order 6, HDF5 output</description>
    <executable>NekMesh</executable>
    <parameters>-m jac:list CubeTet_4_0.msh CubeTet_4_0.nekg:nekg:test</parameters>
    <files>
        <file description="Input File">CubeTet_4_0.msh</file>
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
